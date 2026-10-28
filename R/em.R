# Tempered EM for the beta-binomial mixture at a fixed dispersion nu.
#
# Target, for fixed K:  prod_i sum_k pi_k f(x_i | p_k, nu)^zeta  x  prod_kq Beta(p_kq; a, b),
# f the beta-binomial likelihood of unit i's whole row (binomial coefficient included).
# zeta = 1 is the ordinary mixture likelihood; zeta < 1 flattens each unit's
# contribution (the coarsened posterior of Miller & Dunson 2019), so huge vote
# counts do not split the electorate into ever more blocs. The Beta prior is
# not tempered: the M-step weights are zeta * responsibilities.

#' Fit the BBMM by tempered EM
#'
#' Fits a \eqn{K}-component beta-binomial mixture with a common, fixed
#' dispersion \eqn{\nu} by expectation-maximisation, with each unit's
#' likelihood raised to the power \eqn{\zeta}, from several random starts.
#' The best converged start (highest tempered log-likelihood) is returned.
#'
#' The model: unit \eqn{i} belongs to bloc \eqn{k} with probability
#' \eqn{\pi_k}; given that, its Yes count on question \eqn{q} is
#' \eqn{\mathrm{BetaBinomial}(M_{iq}, \nu p_{kq}, \nu(1 - p_{kq}))}, so
#' \eqn{p_{kq}} is the bloc's expected Yes share and \eqn{\nu} sets how
#' tightly its units cluster around it. Each \eqn{p_{kq}} has a
#' \eqn{\mathrm{Beta}(a, b)} prior (default \eqn{a = b = 10}), so EM finds a
#' posterior mode in \eqn{p} and a maximum-likelihood \eqn{\pi}.
#'
#' Starts: random exponential weights, and the pooled Yes share of each
#' question plus Normal(0, `jitter`) noise for every \eqn{p_{kq}}. With many
#' units EM finds many local optima; `n_starts` of 10 or more is advisable,
#' and `n_at_best` in the result says how many starts reached the best
#' optimum (to within 0.01 in log-likelihood). If only one or two do, the
#' partition may not be unique.
#'
#' @param data A [bbmm_data] object.
#' @param K Number of mixture components.
#' @param zeta Likelihood power in (0, 1]. 1 is the ordinary likelihood.
#' @param nu Fixed beta-binomial dispersion (larger = less over-dispersion).
#' @param prior Shape parameters \eqn{(a, b)} of the Beta prior on each \eqn{p_{kq}}.
#' @param n_starts Number of random starts.
#' @param seed Optional integer seed. The starts are reproducible from it, and
#'   the global random-number state is left untouched. If `NULL`, start seeds
#'   are drawn from the current R random-number stream (so [set.seed()] also
#'   makes runs reproducible).
#' @param tol Convergence tolerance on the change in tempered log-likelihood.
#' @param max_iter Maximum EM iterations per start.
#' @param jitter Standard deviation of the noise in the starting profiles.
#' @param cores Number of processes for the starts ([parallel::mclapply()];
#'   always 1 on Windows).
#' @return An object of class `bbmm_em`, a list with
#' \describe{
#'   \item{K, zeta, nu, prior}{The settings.}
#'   \item{weights}{Mixture weights \eqn{\pi}.}
#'   \item{p}{\eqn{K \times Q} matrix of bloc Yes probabilities.}
#'   \item{gamma}{\eqn{N \times K} matrix of membership probabilities (at power \eqn{\zeta}).}
#'   \item{cluster}{Most probable component of each unit.}
#'   \item{loglik}{Tempered log-likelihood \eqn{\sum_i \log \sum_k \pi_k f_{ik}^\zeta}.}
#'   \item{occupied, min_share}{Components holding at least one unit's
#'     most-probable assignment, and the smallest such share of units.}
#'   \item{BIC}{\eqn{-2\,\mathrm{loglik} + \log(NQ)\,[(T - 1) + TQ]}, \eqn{T} = `occupied`.}
#'   \item{converged, iterations}{For the returned start.}
#'   \item{starts}{Data frame: seed, loglik, converged, iterations of every start.}
#'   \item{n_at_best}{Starts that reached the best optimum.}
#' }
#' @references Miller, J. W. and Dunson, D. B. (2019). Robust Bayesian
#'   inference via coarsening. *Journal of the American Statistical
#'   Association* 114, 1113--1125.
#' @seealso [bbmm_em_select()] to fit a range of \eqn{K}; [bloc_summary()].
#' @examples
#' small <- maine[1:120, 1:20]
#' fit <- bbmm_em(small, K = 3, zeta = 0.3, n_starts = 3, seed = 1)
#' fit
#' @export
bbmm_em <- function(data, K, zeta = 1, nu = 100, prior = c(10, 10), n_starts = 10L, seed = NULL,
                    tol = 1e-5, max_iter = 5000L, jitter = 0.05, cores = 1L) {
  check_data(data)
  K <- as.integer(K)
  if (length(K) != 1 || K < 1) stop("`K` must be one positive integer; see bbmm_em_select() for a range.", call. = FALSE)
  fits <- em_run(data, K, zeta, nu, prior, n_starts, seed, tol, max_iter, jitter, cores)
  fits[[1]]
}

#' Fit the BBMM by tempered EM over a range of K and choose the number of blocs
#'
#' Runs [bbmm_em()] for every \eqn{K} in `K` (all starts in one parallel
#' pool) and applies two rules for the number of blocs:
#' \describe{
#'   \item{BIC}{For each number of occupied components \eqn{T}, the best fit
#'     with that many; the \eqn{T} that minimises BIC.}
#'   \item{Minimum share}{The largest \eqn{K} whose best fit has every occupied
#'     bloc holding at least `min_share` of the units (each unit counted in its
#'     most probable bloc), reported as its number of occupied blocs. With
#'     \eqn{\zeta < 1} this rule is the authors' usual choice (3\% or 5\%).}
#' }
#' BIC with tempered likelihoods is a heuristic: at \eqn{\zeta = 1} with large
#' vote counts it typically keeps rising to the largest \eqn{K} tried.
#'
#' @inheritParams bbmm_em
#' @param K Integer vector of component counts to fit.
#' @param min_share Shares for the minimum-share rule.
#' @return An object of class `bbmm_em_select`: `table` (one row per K: fit
#'   summary, BIC, n_at_best), `selected` (named vector: `BIC`, and one entry
#'   per `min_share`), and `fits` (the `bbmm_em` fit at each K).
#' @seealso [plot_em_selection()]
#' @examples
#' small <- maine[1:120, 1:20]
#' sel <- bbmm_em_select(small, K = 1:4, zeta = 0.3, n_starts = 2, seed = 1)
#' sel
#' @export
bbmm_em_select <- function(data, K = 1:15, zeta = 1, nu = 100, prior = c(10, 10), n_starts = 10L,
                           min_share = c(0.03, 0.05), seed = NULL, tol = 1e-5, max_iter = 5000L,
                           jitter = 0.05, cores = 1L) {
  check_data(data)
  K <- sort(unique(as.integer(K)))
  fits <- em_run(data, K, zeta, nu, prior, n_starts, seed, tol, max_iter, jitter, cores)
  names(fits) <- K
  N <- nrow(data$X); Q <- ncol(data$X)
  tab <- do.call(rbind, lapply(fits, function(f)
    data.frame(K = f$K, loglik = f$loglik, occupied = f$occupied, min_share = f$min_share,
               BIC = f$BIC, converged = f$converged, n_at_best = f$n_at_best)))
  rownames(tab) <- NULL
  by_occ <- tab[order(tab$occupied, -tab$loglik), ]
  by_occ <- by_occ[!duplicated(by_occ$occupied), ]
  sel <- c(BIC = by_occ$occupied[which.min(by_occ$BIC)])
  for (s in min_share) {
    ok <- tab[tab$min_share >= s, ]
    sel[sprintf("min_share_%g%%", 100 * s)] <- if (nrow(ok)) ok$occupied[which.max(ok$K)] else 1L
  }
  structure(list(table = tab, selected = sel, fits = fits, zeta = zeta, nu = nu, N = N, Q = Q,
                 BIC_at_edge = sel[["BIC"]] == max(tab$occupied)), class = "bbmm_em_select")
}

# Internal: all (K, start) EM fits in one pool; returns the best fit per K, in order of K.
em_run <- function(data, Ks, zeta, nu, prior, n_starts, seed, tol, max_iter, jitter, cores) {
  if (!(is.numeric(zeta) && length(zeta) == 1 && zeta > 0 && zeta <= 1)) stop("`zeta` must lie in (0, 1].", call. = FALSE)
  if (!(nu > 0)) stop("`nu` must be positive.", call. = FALSE)
  if (length(prior) != 2 || any(prior <= 0)) stop("`prior` must be two positive numbers.", call. = FALSE)
  n_starts <- as.integer(n_starts)
  x <- data$X; n <- data$M
  miss <- is.na(x) | is.na(n); x[miss] <- 0L; n[miss] <- 0L
  storage.mode(x) <- storage.mode(n) <- "double"
  seeds <- start_seeds(seed, length(Ks) * n_starts)
  jobs <- data.frame(K = rep(Ks, each = n_starts), start = rep(seq_len(n_starts), length(Ks)), seed = seeds)
  jobs <- jobs[order(-jobs$K, jobs$start), ]                       # big K first: better load balance
  one <- function(j) em_one(x, n, jobs$K[j], zeta, nu, prior, jobs$seed[j], tol, max_iter, jitter)
  res <- run_parallel(seq_len(nrow(jobs)), one, cores)
  lapply(Ks, function(K) {
    r <- res[jobs$K == K]
    starts <- data.frame(seed = vapply(r, `[[`, 1, "seed"), loglik = vapply(r, `[[`, 1, "loglik"),
                         converged = vapply(r, `[[`, TRUE, "converged"),
                         iterations = vapply(r, `[[`, 1L, "iterations"))
    pool <- if (any(starts$converged)) which(starts$converged) else seq_along(r)
    b <- pool[which.max(starts$loglik[pool])]
    f <- r[[b]]
    f$starts <- starts[order(-starts$loglik), ]
    f$n_at_best <- sum(starts$loglik[pool] > f$loglik - 0.01)
    em_finish(f, data, zeta, nu, prior)
  })
}

# Internal: one EM run from one random start.
em_one <- function(x, n, K, zeta, nu, prior, seed, tol, max_iter, jitter) {
  N <- nrow(x); Q <- ncol(x)
  with_seed(seed, {
    w <- rexp(K); w <- w / sum(w)
    pooled <- colSums(x) / pmax(colSums(n), 1)
    p <- matrix(rep(pooled, each = K), K, Q) + matrix(rnorm(K * Q, sd = jitter), K, Q)
  })
  p <- pmin(pmax(p, 1e-3), 1 - 1e-3)
  L <- zeta * bb_loglik_em_cpp(x, n, p, nu)
  e <- e_step(L, w); ll <- e$ll; converged <- FALSE; it <- 0L
  while (it < max_iter) {
    it <- it + 1L
    w <- colSums(e$gamma) / N
    p <- bb_newton_em_cpp(zeta * e$gamma, x, n, nu, p, prior[1], prior[2], 75L, 1e-8, 1e-6)
    L <- zeta * bb_loglik_em_cpp(x, n, p, nu)
    e <- e_step(L, w)
    if (!is.finite(e$ll)) break
    if (abs(e$ll - ll) < tol) { ll <- e$ll; converged <- TRUE; break }
    ll <- e$ll
  }
  list(K = K, seed = seed, weights = w, p = p, gamma = e$gamma, loglik = ll, converged = converged, iterations = it)
}

e_step <- function(L, w) {
  lj <- sweep(L, 2L, log(pmax(w, 1e-300)), "+")
  lse <- lse_rows(lj)
  list(gamma = exp(lj - lse), ll = sum(lse))
}

lse_rows <- function(A) {
  m <- apply(A, 1L, max)
  m + log(rowSums(exp(A - m)))
}

em_finish <- function(f, data, zeta, nu, prior) {
  N <- nrow(data$X); Q <- ncol(data$X)
  dimnames(f$p) <- list(paste0("component_", seq_len(f$K)), data$questions)
  dimnames(f$gamma) <- list(data$units, rownames(f$p))
  f$cluster <- max.col(f$gamma, ties.method = "first")
  sh <- tabulate(f$cluster, f$K) / N
  f$occupied <- sum(sh > 0)
  f$min_share <- min(sh[sh > 0])
  f$BIC <- -2 * f$loglik + log(N * Q) * ((f$occupied - 1) + f$occupied * Q)
  f$zeta <- zeta; f$nu <- nu; f$prior <- prior
  f$units <- data$units; f$questions <- data$questions
  class(f) <- "bbmm_em"
  f
}

#' @export
print.bbmm_em <- function(x, ...) {
  cat(sprintf("BBMM tempered EM: K = %d (%d occupied), zeta = %g, nu = %g\n", x$K, x$occupied, x$zeta, x$nu))
  cat(sprintf("Tempered log-likelihood %.2f, BIC %.1f; %s after %d iterations; %d of %d starts reached the best optimum\n",
              x$loglik, x$BIC, if (x$converged) "converged" else "NOT converged", x$iterations,
              x$n_at_best, nrow(x$starts)))
  sizes <- tabulate(x$cluster, x$K)
  cat("Units per component:", paste(sizes, collapse = " "), "\n")
  invisible(x)
}

#' @export
print.bbmm_em_select <- function(x, ...) {
  cat(sprintf("BBMM tempered EM over K = %s (zeta = %g, nu = %g; %d units x %d questions)\n",
              paste(range(x$table$K), collapse = "-"), x$zeta, x$nu, x$N, x$Q))
  print(x$table, row.names = FALSE, digits = 6)
  cat("\nSelected number of blocs:\n")
  print(x$selected)
  if (isTRUE(x$BIC_at_edge)) cat("Note: BIC is still falling at the largest K tried.\n")
  invisible(x)
}

# ---- helpers shared by the engines ------------------------------------------

# Seeds for n jobs: from `seed` if given (without touching the global stream), else from R's stream.
start_seeds <- function(seed, n) {
  if (is.null(seed)) return(sample.int(.Machine$integer.max - 1L, n))
  with_seed(seed, sample.int(.Machine$integer.max - 1L, n))
}

# Evaluate expr with set.seed(seed), restoring the caller's random-number state afterwards.
with_seed <- function(seed, expr) {
  env <- globalenv()
  old <- if (exists(".Random.seed", envir = env, inherits = FALSE)) get(".Random.seed", envir = env) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = env) else assign(".Random.seed", old, envir = env))
  set.seed(seed)
  expr
}

run_parallel <- function(idx, fun, cores) {
  cores <- as.integer(cores)
  if (.Platform$OS.type == "windows" || cores <= 1L) return(lapply(idx, fun))
  parallel::mclapply(idx, fun, mc.cores = cores, mc.preschedule = FALSE)
}
