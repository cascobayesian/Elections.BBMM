# Fixed-K MCMC at likelihood power zeta: the stepping-stone kernel run as chains.

#' Fit the Bayesian BBMM at fixed K by tempered MCMC
#'
#' Samples the posterior of the \eqn{K}-component Bayesian BBMM with the
#' likelihood raised to the power \eqn{\zeta} (the coarsened posterior), by
#' running the stepping-stone kernel ([bbmm_ss()]) as ordinary Markov chains
#' at \eqn{\beta = \zeta}. Each iteration is `n_sweeps` sweeps of: redraw empty
#' components from the prior; `n_sm` Jain--Neal split-merge proposals; a
#' collapsed Gibbs sweep of the allocations (\eqn{\pi} integrated out);
#' random-walk Metropolis on every logit \eqn{p_{kq}} and log \eqn{\nu_k}.
#' Unlike [bbmm_em()], every bloc has its own dispersion \eqn{\nu_k} with a
#' Gamma prior, and components may empty out, so the number of *occupied*
#' blocs \eqn{T \le K} is part of the posterior.
#'
#' Starting from the EM fit (the default) puts every chain near a good mode at
#' once; the chains then explore around it. Chains started from different EM
#' optima (`init = "em"` uses the best `n_chains` distinct starts) that never
#' agree are a sign that the posterior is multimodal and that the summaries
#' describe the neighbourhood of one mode. The random-walk sds adapt during
#' burn-in only, so the saved draws come from a fixed kernel.
#'
#' **Run time.** Roughly the cost of one SMC step per iteration: on Maine
#' (462 x 72), \eqn{K = 6}, about 10--30 iterations per second per chain.
#'
#' @inheritParams bbmm_ss
#' @param n_iter Iterations saved per chain (after burn-in), before thinning.
#' @param n_burn Burn-in iterations per chain.
#' @param thin Keep every `thin`-th iteration.
#' @param n_chains Number of chains (run in parallel over `n_threads`).
#' @param init `"em"` (fit [bbmm_em()] at this K and zeta and start each chain
#'   from one of its best starts), `"prior"` (a draw from the prior), a
#'   `bbmm_em` object, or a vector of starting allocations in `1:K`.
#' @param control Settings, see [bbmm_control()]; `n_sweeps` defaults to 1 here.
#' @return An object of class `bbmm_mcmc`: draws `z` (draws x units), `p`
#'   (draws x K x Q), `nu` (draws x K), `ll` (untempered complete-data
#'   log-likelihood), `T` (occupied components), `chain` and `iteration` of each
#'   draw, `acceptance` (per chain, after burn-in) and the settings.
#' @seealso [bloc_summary()], [plot_trace()], [mcmc_diagnostics()]
#' @examples
#' tiny <- maine[1:60, 1:10]
#' fit <- bbmm_mcmc(tiny, K = 3, zeta = 0.3, n_iter = 100, n_burn = 50, n_chains = 2,
#'                  seed = 1, progress = FALSE)
#' fit
#' @export
bbmm_mcmc <- function(data, K, zeta = 1, n_iter = 2000L, n_burn = 1000L, thin = 1L, n_chains = 2L,
                      init = "em", control = list(), seed = NULL, n_threads = 1L, progress = interactive()) {
  check_data(data)
  K <- as.integer(K); n_chains <- as.integer(n_chains)
  if (is.null(control$n_sweeps)) control$n_sweeps <- 1L
  seed <- cpp_seed(seed)
  sc <- cpp_control(K, control, zeta, seed, n_threads, progress, n_particles = 2L)
  N <- nrow(data$X); Q <- ncol(data$X)
  st <- mcmc_start(data, K, zeta, init, n_chains, sc$ctl, seed, n_threads)
  cm <- cpp_counts(data)
  t0 <- proc.time()[["elapsed"]]
  out <- bbmm_mcmc_cpp(cm$X, cm$M, sc$cc, zeta, as.integer(n_iter), as.integer(n_burn), as.integer(thin),
                       st$z, st$p, st$nu, st$mu, isTRUE(progress))
  out$seconds <- proc.time()[["elapsed"]] - t0
  dimnames(out$p) <- list(NULL, NULL, data$questions)
  out$T <- apply(out$z, 1, function(z) length(unique(z)))
  if (!isTRUE(sc$ctl$nu_pool)) out$mu <- NULL
  rownames(out$acceptance) <- paste0("chain_", seq_len(n_chains))
  out$weights <- rep(1 / nrow(out$z), nrow(out$z))
  out[c("K", "zeta", "n_iter", "n_burn", "thin", "n_chains", "control", "seed", "init")] <-
    list(K, zeta, n_iter, n_burn, thin, n_chains, sc$ctl, seed, st$label)
  out$units <- data$units; out$questions <- data$questions
  class(out) <- "bbmm_mcmc"
  out
}

# Internal: starting states, n_chains x N allocations, n_chains x K x Q profiles, n_chains x K dispersions.
mcmc_start <- function(data, K, zeta, init, n_chains, ctl, seed, n_threads) {
  N <- nrow(data$X); Q <- ncol(data$X)
  nu0 <- max(100, ctl$nu_min)
  pooled_p <- function(z) {
    X <- data$X; M <- data$M; X[is.na(X)] <- 0L; M[is.na(M)] <- 0L
    p <- matrix(0.5, K, Q)
    for (k in unique(z)) p[k, ] <- (colSums(X[z == k, , drop = FALSE]) + ctl$p_prior_a) /
      (colSums(M[z == k, , drop = FALSE]) + ctl$p_prior_a + ctl$p_prior_b)
    p
  }
  Z <- matrix(0L, n_chains, N); P <- array(0, c(n_chains, K, Q)); NU <- matrix(nu0, n_chains, K)
  if (identical(init, "em") || inherits(init, "bbmm_em")) {
    if (identical(init, "em")) {
      init <- em_run(data, K, zeta, 100, c(ctl$p_prior_a, ctl$p_prior_b), max(5L, n_chains), seed,
                     1e-5, 5000L, 0.05, 1L)[[1]]
      label <- "em"
    } else label <- "bbmm_em object"
    if (init$K != K) stop("The bbmm_em fit has K = ", init$K, ", not ", K, ".", call. = FALSE)
    # chain c starts from the c-th best start: rerun those starts to recover their states
    st <- init$starts
    use <- head(st$seed[!duplicated(round(st$loglik, 2))], n_chains)
    use <- rep_len(use, n_chains)
    x <- data$X; n <- data$M; miss <- is.na(x) | is.na(n); x[miss] <- 0L; n[miss] <- 0L
    storage.mode(x) <- storage.mode(n) <- "double"
    for (c in seq_len(n_chains)) {
      f <- if (use[c] == init$seed) init else
        em_one(x, n, K, init$zeta, init$nu, init$prior, use[c], 1e-5, 5000L, 0.05)
      Z[c, ] <- max.col(f$gamma, ties.method = "first")
      P[c, , ] <- f$p
    }
  } else if (identical(init, "prior")) {
    label <- "prior"
    with_seed(seed, {
      for (c in seq_len(n_chains)) {
        Z[c, ] <- sample.int(K, N, replace = TRUE)
        P[c, , ] <- rbeta(K * Q, ctl$p_prior_a, ctl$p_prior_b)
        nus <- rgamma(10 * K, ctl$nu_shape, rate = ctl$nu_rate)
        NU[c, ] <- head(c(nus[nus >= ctl$nu_min], rep(nu0, K)), K)
      }
    })
  } else if (is.numeric(init) && length(init) == N) {
    label <- "given allocations"
    z <- as.integer(init)
    if (any(z < 1 | z > K)) stop("Starting allocations must lie in 1..K.", call. = FALSE)
    for (c in seq_len(n_chains)) { Z[c, ] <- z; P[c, , ] <- pooled_p(z) }
  } else stop("`init` must be \"em\", \"prior\", a bbmm_em fit or a vector of allocations.", call. = FALSE)
  if (ctl$nu_pool) NU[] <- exp(ctl$pool_m0)
  storage.mode(Z) <- "integer"
  list(z = Z, p = P, nu = NU, mu = if (ctl$nu_pool) rep(ctl$pool_m0, n_chains) else numeric(0), label = label)
}

#' @export
print.bbmm_mcmc <- function(x, ...) {
  cat(sprintf("BBMM fixed-K MCMC: K = %d, zeta = %g; %d chain(s) x %d saved draws (burn-in %d, thin %d), %.0f s\n",
              x$K, x$zeta, x$n_chains, nrow(x$z) / x$n_chains, x$n_burn, x$thin, x$seconds))
  tt <- table(x$T) / length(x$T)
  cat("Posterior on occupied components:", paste(sprintf("T=%s: %.2f", names(tt), tt), collapse = ", "), "\n")
  cat("Acceptance after burn-in (mean over chains):",
      paste(sprintf("%s %.2f", colnames(x$acceptance), colMeans(x$acceptance, na.rm = TRUE)), collapse = ", "), "\n")
  d <- mcmc_diagnostics(x)
  if (!is.na(d$rhat_ll)) cat(sprintf("Split R-hat of the log-likelihood: %.3f\n", d$rhat_ll))
  invisible(x)
}

#' Convergence diagnostics for bbmm_mcmc
#'
#' Split-\eqn{\hat R} (Gelman et al. 2013) of the log-likelihood and of the
#' number of occupied components across chains, each chain's mean
#' log-likelihood, and the adjusted Rand index between every pair of chains'
#' final allocations. Values of \eqn{\hat R} well above 1.01 or low
#' between-chain agreement mean the chains sit in different modes; a chain
#' whose mean log-likelihood is far below the others is stuck in a worse one.
#'
#' @param fit A `bbmm_mcmc` object.
#' @return A list: `rhat_ll`, `rhat_T`, `mean_ll` (per chain), `ari` (chains x
#'   chains matrix of adjusted Rand indices between final allocations).
#' @examples
#' tiny <- maine[1:60, 1:10]
#' fit <- bbmm_mcmc(tiny, K = 3, zeta = 0.3, n_iter = 100, n_burn = 50, seed = 1, progress = FALSE)
#' mcmc_diagnostics(fit)
#' @export
mcmc_diagnostics <- function(fit) {
  stopifnot(inherits(fit, "bbmm_mcmc"))
  last <- vapply(split(seq_along(fit$chain), fit$chain), max, 1L)
  ari <- outer(seq_along(last), seq_along(last),
               Vectorize(function(a, b) adjusted_rand(fit$z[last[a], ], fit$z[last[b], ])))
  dimnames(ari) <- list(paste0("chain_", seq_along(last)), paste0("chain_", seq_along(last)))
  list(rhat_ll = split_rhat(fit$ll, fit$chain), rhat_T = split_rhat(fit$T, fit$chain),
       mean_ll = tapply(fit$ll, fit$chain, mean), ari = ari)
}

split_rhat <- function(v, chain) {
  parts <- unlist(lapply(split(v, chain), function(x) {
    h <- floor(length(x) / 2); if (h < 2) return(NULL); list(x[1:h], x[(h + 1):(2 * h)]) }), recursive = FALSE)
  if (length(parts) < 2) return(NA_real_)
  n <- length(parts[[1]]); m <- vapply(parts, mean, 1); s2 <- vapply(parts, var, 1)
  W <- mean(s2); B <- n * var(m)
  if (W == 0) return(if (B == 0) 1 else Inf)
  sqrt(((n - 1) / n * W + B / n) / W)
}
