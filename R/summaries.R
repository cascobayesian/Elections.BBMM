# Posterior summaries of the blocs, common to the three engines.
#
# EM gives one partition and membership probabilities directly. The samplers
# give many partitions whose labels mean nothing across draws (label
# switching), so their summaries start from label-free quantities:
#   - the posterior similarity matrix psm[i, j] = P(i and j share a bloc | x);
#   - a representative partition: the sampled partition closest to psm in
#     squared error (Dahl 2006, "least-squares clustering");
# and only then attach parameters to the representative blocs, by matching
# each draw's labels to them (greatest overlap first).
# Blocs are numbered by their mean Yes share across questions, highest first.

#' Summarise the blocs of a fitted BBMM
#'
#' Turns any fit into the same set of tables: which bloc each unit belongs to
#' and how surely, how big each bloc is, and how each bloc voted on each
#' question.
#'
#' For [bbmm_em()] fits the partition is each unit's most probable component
#' and `certainty` its membership probability. For [bbmm_mcmc()],
#' [bbmm_ss()] and [bbmm_evidence()] fits:
#' \itemize{
#'   \item `psm`, the posterior similarity matrix: the posterior probability
#'     that two units are in the same bloc;
#'   \item the representative partition is the sampled partition that
#'     minimises \eqn{\sum_{ij} (1[z_i = z_j] - psm_{ij})^2}
#'     (least-squares clustering, Dahl 2006), optionally among draws with
#'     exactly `n_blocs` occupied blocs;
#'   \item a unit's `certainty` is its average posterior probability of
#'     sharing a bloc with the other members of its representative bloc, and
#'     `membership[i, b]` the same average for every bloc `b`;
#'   \item bloc parameters (\eqn{p_{kq}}, \eqn{\nu_k}) are summarised by
#'     matching each draw's components to the representative blocs (largest
#'     overlap first) and taking weighted means and 90\% intervals.
#' }
#' For a `bbmm_evidence` object the runs at each \eqn{K} are pooled with
#' weights \eqn{p(K \mid x)} times the particle weights.
#'
#' @param fit A `bbmm_em`, `bbmm_mcmc`, `bbmm_ss` or `bbmm_evidence` object.
#' @param data The [bbmm_data] the fit used.
#' @param n_blocs Optional: choose the representative partition among draws
#'   with exactly this many occupied blocs (samplers only).
#' @param ... Unused.
#' @return An object of class `bloc_summary`, a list with
#' \describe{
#'   \item{assignment}{Data frame: `unit`, `bloc`, `certainty`.}
#'   \item{membership}{Units x blocs matrix (see Details).}
#'   \item{blocs}{Data frame per bloc: `n_units`, `share` of units, `ballots`
#'     (total two-option votes cast), `mean_yes` (mean over questions of the
#'     bloc's ballot-weighted Yes share), `nu` (dispersion; posterior median
#'     and 90\% interval for the samplers), `certainty` (mean over members).}
#'   \item{questions}{Long data frame per bloc and question: `yes_share`
#'     (ballot-weighted, from the data), `all_share` (all units), `gap`
#'     (= yes_share - all_share), and the model's `p_mean`, `p_lo`, `p_hi`.}
#'   \item{psm}{Posterior similarity matrix (samplers only).}
#'   \item{n_blocs}{Posterior distribution of the number of occupied blocs (samplers only).}
#' }
#' @references Dahl, D. B. (2006). Model-based clustering for expression data
#'   via a Dirichlet process mixture model. In *Bayesian Inference for Gene
#'   Expression and Proteomics*, Cambridge University Press, 201--218.
#' @seealso [distinctive_questions()], [question_fit()], [plot_bloc_map()],
#'   [plot_bloc_profiles()], [plot_js_pca()], [plot_coclustering()]
#' @examples
#' small <- maine[1:120, 1:20]
#' fit <- bbmm_em(small, K = 3, zeta = 0.3, n_starts = 3, seed = 1)
#' s <- bloc_summary(fit, small)
#' s
#' head(s$assignment)
#' @export
bloc_summary <- function(fit, data, ...) UseMethod("bloc_summary")

#' @rdname bloc_summary
#' @export
bloc_summary.bbmm_em <- function(fit, data, ...) {
  check_match(fit, data)
  occ <- sort(unique(fit$cluster))
  z <- match(fit$cluster, occ)
  ord <- bloc_order(z, data)
  z <- ord[z]
  comp <- occ[order(ord)]                                  # component behind bloc b
  membership <- fit$gamma[, comp, drop = FALSE]
  colnames(membership) <- paste0("bloc_", seq_along(comp))
  B <- length(comp); Q <- ncol(data$X)
  pm <- fit$p[comp, , drop = FALSE]
  out <- assemble_summary(z, data, certainty = membership[cbind(seq_along(z), z)], membership = membership,
                          p_mean = pm, p_lo = pm, p_hi = pm,
                          nu = matrix(fit$nu, B, 3, dimnames = list(NULL, c("median", "lo", "hi"))))
  out$method <- sprintf("tempered EM, K = %d (%d occupied), zeta = %g, fixed nu = %g", fit$K, B, fit$zeta, fit$nu)
  out
}

#' @rdname bloc_summary
#' @export
bloc_summary.bbmm_mcmc <- function(fit, data, n_blocs = NULL, ...) {
  check_match(fit, data)
  out <- draws_summary(fit$z, fit$weights, list(p = fit$p, nu = fit$nu), rep(1L, nrow(fit$z)), data, n_blocs)
  out$method <- sprintf("fixed-K MCMC, K = %d, zeta = %g, %d draws from %d chain(s)", fit$K, fit$zeta,
                        nrow(fit$z), fit$n_chains)
  out
}

#' @rdname bloc_summary
#' @export
bloc_summary.bbmm_ss <- function(fit, data, n_blocs = NULL, ...) {
  check_match(fit, data)
  if (is.null(fit$z)) stop("The run kept no particles (keep_particles = FALSE).", call. = FALSE)
  out <- draws_summary(fit$z, fit$weights, list(p = fit$p, nu = fit$nu), rep(1L, nrow(fit$z)), data, n_blocs)
  out$method <- sprintf("stepping-stone SMC, K = %d, zeta = %g, %d particles", fit$K, fit$zeta, nrow(fit$z))
  out
}

#' @rdname bloc_summary
#' @export
bloc_summary.bbmm_evidence <- function(fit, data, n_blocs = NULL, ...) {
  if (is.null(fit$fits)) stop("The evidence object kept no fits (keep_fits = FALSE).", call. = FALSE)
  check_match(fit$fits[[1]], data)
  pk <- fit$posterior_K
  Z <- list(); W <- list(); src <- list(); pars <- list()
  for (s in seq_along(fit$fits)) {
    f <- fit$fits[[s]]
    wk <- pk$posterior[pk$K == f$K] / pk$n_rep[pk$K == f$K]
    if (wk < 1e-6) next
    Z[[s]] <- f$z; W[[s]] <- wk * f$weights; src[[s]] <- rep(s, nrow(f$z))
    pars[[s]] <- list(p = f$p, nu = f$nu)
  }
  keep <- !vapply(Z, is.null, TRUE)
  out <- draws_summary(do.call(rbind, Z[keep]), unlist(W[keep]), pars, unlist(src[keep]), data, n_blocs)
  Ks <- sort(unique(fit$runs$K))
  k_lab <- if (length(Ks) > 2 && all(diff(Ks) == 1)) paste0(Ks[1], "-", Ks[length(Ks)]) else paste(Ks, collapse = ",")
  out$method <- sprintf("stepping-stone SMC pooled over K = %s by p(K | x), zeta = %g", k_lab, fit$zeta)
  out
}

# Internal: summaries from weighted draws of partitions.
#   Z: draws x N; w: weights; pars: list (indexed by src) of list(p = draws x K x Q, nu = draws x K)
#   with rows in the order the draws of that source appear in Z; src: source of each draw.
draws_summary <- function(Z, w, pars, src, data, n_blocs) {
  N <- ncol(Z)
  if (!is.null(pars$p)) pars <- list(pars)
  row_in_src <- stats::ave(seq_along(src), src, FUN = seq_along)
  w <- w / sum(w)
  Tz <- apply(Z, 1, function(z) length(unique(z)))
  n_blocs_post <- tapply(w, Tz, sum)
  # identical partitions (common after SMC resampling) are scored once
  key <- apply(Z, 1, function(z) paste(match(z, unique(z)), collapse = ","))
  uk <- !duplicated(key)
  wu <- tapply(w, factor(key, levels = key[uk]), sum)
  Zu <- Z[uk, , drop = FALSE]
  psm <- matrix(0, N, N)
  for (r in seq_len(nrow(Zu))) {
    A <- onehot(Zu[r, ])
    psm <- psm + wu[[r]] * tcrossprod(A)
  }
  dimnames(psm) <- list(data$units, data$units)
  cand <- seq_len(nrow(Zu))
  if (!is.null(n_blocs)) {
    cand <- which(Tz[uk] == n_blocs)
    if (!length(cand)) stop("No draw has exactly ", n_blocs, " occupied blocs.", call. = FALSE)
  }
  loss <- vapply(cand, function(r) {
    A <- onehot(Zu[r, ]); nk <- colSums(A)
    sum(nk^2) - 2 * sum(A * (psm %*% A))
  }, 1)
  zr <- Zu[cand[which.min(loss)], ]
  z <- match(zr, sort(unique(zr)))
  ord <- bloc_order(z, data); z <- ord[z]
  B <- max(z)
  A <- onehot(z, B)
  nb <- colSums(A)
  S <- psm %*% A - A * diag(psm)                         # summed similarity to each bloc's other members
  denom <- matrix(nb, N, B, byrow = TRUE) - A
  membership <- ifelse(denom > 0, S / pmax(denom, 1), NA_real_)
  colnames(membership) <- paste0("bloc_", seq_len(B)); rownames(membership) <- data$units
  certainty <- membership[cbind(seq_len(N), z)]
  certainty[is.na(certainty)] <- 1
  # parameters: match each draw's components to the representative blocs
  Q <- ncol(data$X)
  pd <- array(NA_real_, c(nrow(Z), B, Q)); nud <- matrix(NA_real_, nrow(Z), B)
  for (r in seq_len(nrow(Z))) {
    m <- match_labels(Z[r, ], z, B)
    pr <- pars[[src[r]]]
    ok <- !is.na(m)
    pd[r, ok, ] <- pr$p[row_in_src[r], m[ok], ]
    nud[r, ok] <- pr$nu[row_in_src[r], m[ok]]
  }
  summ <- function(v) { ok <- !is.na(v); if (!any(ok)) return(c(NA, NA, NA, NA))
    c(sum(w[ok] * v[ok]) / sum(w[ok]), wquantile(v[ok], w[ok], c(0.05, 0.5, 0.95))) }
  ps <- apply(pd, c(2, 3), summ)                        # 4 x B x Q
  nus <- t(apply(nud, 2, summ))
  out <- assemble_summary(z, data, certainty, membership, p_mean = ps[1, , ], p_lo = ps[2, , ], p_hi = ps[4, , ],
                          nu = cbind(median = nus[, 3], lo = nus[, 2], hi = nus[, 4]))
  out$psm <- psm
  out$n_blocs <- data.frame(T = as.integer(names(n_blocs_post)), posterior = as.numeric(n_blocs_post))
  out$n_draws <- nrow(Z); out$n_distinct <- nrow(Zu)
  out
}

# Internal: build the tables from a partition z (1..B, already ordered).
assemble_summary <- function(z, data, certainty, membership, p_mean, p_lo, p_hi, nu) {
  B <- max(z); N <- length(z)
  X <- data$X; M <- data$M; X[is.na(X)] <- 0L; M[is.na(M)] <- 0L
  Xb <- rowsum(X, z, reorder = TRUE); Mb <- rowsum(M, z, reorder = TRUE)
  share <- Xb / Mb
  all <- colSums(X) / colSums(M)
  p_mean <- matrix(p_mean, B, ncol(X)); p_lo <- matrix(p_lo, B, ncol(X)); p_hi <- matrix(p_hi, B, ncol(X))
  qs <- data.frame(bloc = rep(seq_len(B), ncol(X)), question = rep(data$questions, each = B),
                   label = rep(data$labels, each = B), yes_share = as.vector(share), all_share = rep(all, each = B),
                   gap = as.vector(share - matrix(all, B, ncol(X), byrow = TRUE)),
                   p_mean = as.vector(p_mean), p_lo = as.vector(p_lo), p_hi = as.vector(p_hi))
  qs <- qs[order(qs$bloc, match(qs$question, data$questions)), ]
  rownames(qs) <- NULL
  blocs <- data.frame(bloc = seq_len(B), n_units = tabulate(z, B), share = tabulate(z, B) / N,
                      ballots = rowSums(Mb), mean_yes = rowMeans(share, na.rm = TRUE),
                      nu = nu[, 1], nu_lo = nu[, 2], nu_hi = nu[, 3],
                      certainty = as.numeric(tapply(certainty, factor(z, levels = seq_len(B)), mean)))
  structure(list(assignment = data.frame(unit = data$units, bloc = z, certainty = certainty),
                 membership = membership, blocs = blocs, questions = qs, psm = NULL, n_blocs = NULL),
            class = "bloc_summary")
}

# Internal: map labels 1..B to new numbers ordered by mean Yes share (highest = 1).
bloc_order <- function(z, data) {
  X <- data$X; M <- data$M; X[is.na(X)] <- 0L; M[is.na(M)] <- 0L
  sh <- rowsum(X, z, reorder = TRUE) / rowsum(M, z, reorder = TRUE)
  m <- rowMeans(sh, na.rm = TRUE)
  rank(-m, ties.method = "first")
}

onehot <- function(z, B = max(z)) {
  A <- matrix(0, length(z), B)
  A[cbind(seq_along(z), z)] <- 1
  A
}

# Internal: for each representative bloc b, the draw's label with the largest
# overlap, matched greedily one-to-one (NA when the draw has too few labels).
match_labels <- function(zd, zr, B) {
  O <- table(factor(zr, levels = seq_len(B)), zd)
  labs <- as.integer(colnames(O)); O <- unclass(O)
  m <- rep(NA_integer_, B)
  while (any(O > 0)) {
    ij <- which(O == max(O), arr.ind = TRUE)[1, ]
    m[ij[1]] <- labs[ij[2]]
    O[ij[1], ] <- -1; O[, ij[2]] <- -1
  }
  m
}

wquantile <- function(v, w, probs) {
  o <- order(v); v <- v[o]; cw <- cumsum(w[o]) / sum(w)
  vapply(probs, function(p) v[min(which(cw >= p - 1e-12))], 1)
}

check_match <- function(fit, data) {
  check_data(data)
  if (!identical(fit$units, data$units) || !identical(fit$questions, data$questions))
    stop("`data` does not match the data the fit used (units or questions differ).", call. = FALSE)
}

#' Adjusted Rand index between two partitions
#'
#' @param a,b Vectors of labels of the same length.
#' @return The adjusted Rand index (1 = identical partitions, about 0 = chance agreement).
#' @examples
#' adjusted_rand(c(1, 1, 2, 2, 3), c(2, 2, 1, 1, 1))
#' @export
adjusted_rand <- function(a, b) {
  tab <- table(a, b); c2 <- function(n) n * (n - 1) / 2
  s <- sum(c2(tab)); sa <- sum(c2(rowSums(tab))); sb <- sum(c2(colSums(tab)))
  e <- sa * sb / c2(length(a))
  if ((sa + sb) / 2 == e) return(1)
  (s - e) / ((sa + sb) / 2 - e)
}

#' @export
print.bloc_summary <- function(x, ...) {
  cat("Bloc summary:", x$method, "\n")
  b <- x$blocs
  print(data.frame(bloc = b$bloc, units = b$n_units, share = round(b$share, 3),
                   ballots = format(b$ballots, big.mark = ","), mean_yes = round(b$mean_yes, 3),
                   nu = round(b$nu), certainty = round(b$certainty, 3)), row.names = FALSE)
  if (!is.null(x$n_blocs)) {
    nb <- x$n_blocs[x$n_blocs$posterior > 0.005, ]
    cat("Posterior on the number of occupied blocs:", paste(sprintf("%d: %.2f", nb$T, nb$posterior), collapse = ", "), "\n")
  }
  invisible(x)
}

#' The questions that set each bloc apart
#'
#' For each bloc, the questions on which its ballot-weighted Yes share differs
#' most from that of all units together.
#'
#' @param summary A [bloc_summary] object.
#' @param n Questions per bloc.
#' @return Data frame: `bloc`, `label`, `yes_share`, `all_share`, `gap`, ordered by bloc and |gap|.
#' @examples
#' small <- maine[1:120, 1:20]
#' s <- bloc_summary(bbmm_em(small, K = 3, zeta = 0.3, n_starts = 3, seed = 1), small)
#' distinctive_questions(s, n = 3)
#' @export
distinctive_questions <- function(summary, n = 5L) {
  stopifnot(inherits(summary, "bloc_summary"))
  q <- summary$questions
  out <- do.call(rbind, lapply(split(q, q$bloc), function(t) head(t[order(-abs(t$gap)), ], n)))
  rownames(out) <- NULL
  out[c("bloc", "question", "label", "yes_share", "all_share", "gap")]
}

#' How well the blocs account for each question
#'
#' For each question: the spread between the blocs' Yes shares, and the
#' share of the variation in units' Yes shares that bloc membership explains
#' (ballot-weighted \eqn{R^2}). Questions with a low \eqn{R^2} cut across the
#' blocs -- the units' votes on them follow some other pattern -- and are the
#' ones to read closely.
#'
#' @param summary A [bloc_summary] object.
#' @param data The [bbmm_data] used for the fit.
#' @return Data frame per question: `label`, `all_share`, `spread_data`
#'   (max - min of the blocs' ballot-weighted Yes shares), `spread_model` (max
#'   - min of the blocs' model Yes probabilities), `r2` (ballot-weighted),
#'   ordered by `r2`.
#' @examples
#' small <- maine[1:120, 1:20]
#' s <- bloc_summary(bbmm_em(small, K = 3, zeta = 0.3, n_starts = 3, seed = 1), small)
#' head(question_fit(s, small))
#' @export
question_fit <- function(summary, data) {
  stopifnot(inherits(summary, "bloc_summary"))
  check_data(data)
  if (!identical(summary$assignment$unit, data$units)) stop("`data` does not match the summary.", call. = FALSE)
  z <- summary$assignment$bloc
  q <- summary$questions
  out <- do.call(rbind, lapply(seq_along(data$questions), function(j) {
    x <- data$X[, j]; m <- data$M[, j]; ok <- !is.na(m)
    y <- x[ok] / m[ok]; wt <- m[ok]
    t <- q[q$question == data$questions[j], ]
    yb <- t$yes_share[match(z[ok], t$bloc)]
    yall <- sum(x[ok]) / sum(wt)
    data.frame(question = data$questions[j], label = data$labels[j], all_share = yall,
               spread_data = diff(range(t$yes_share, na.rm = TRUE)),
               spread_model = diff(range(t$p_mean, na.rm = TRUE)),
               r2 = 1 - sum(wt * (y - yb)^2) / sum(wt * (y - yall)^2))
  }))
  out <- out[order(out$r2), ]
  rownames(out) <- NULL
  out
}
