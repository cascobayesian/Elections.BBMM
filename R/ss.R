# Bayesian BBMM at fixed K: controls shared by the stepping-stone SMC (bbmm_ss)
# and the fixed-K MCMC (bbmm_mcmc), both of which run the kernel in src/ss_smc.cpp.
#
# Model (K components; empty components allowed):
#   pi ~ Dirichlet(gamma, ..., gamma), integrated out
#   p_kq ~ Beta(p_prior_a, p_prior_b)
#   nu_k ~ Gamma(nu_shape, rate nu_rate) truncated to [nu_min, Inf), one per component
#   X_iq | z_i = k ~ BetaBinomial(M_iq, nu_k p_kq, nu_k (1 - p_kq)); missing cells contribute nothing
# Target at likelihood power zeta: prior x L(x | z, p, nu)^zeta. The binomial
# coefficient is omitted, so log Z differs from the full one by a data-only
# constant that is the same for every K.

#' Control settings for the Bayesian samplers
#'
#' Settings shared by [bbmm_ss()] and [bbmm_mcmc()]. Call with no arguments
#' to see every default. The defaults are the settings the authors use for
#' production runs (200 particles, conditional ESS 0.9, two kernel sweeps and
#' five Jain--Neal split-merge proposals per SMC step, \eqn{\nu \ge 50}).
#'
#' @param n_particles Number of SMC particles.
#' @param cess_target Conditional ESS fraction that sets each tempering step
#'   (larger = smaller steps, slower, more accurate).
#' @param resample_threshold Resample when ESS < threshold x particles.
#' @param n_sweeps Kernel sweeps per particle per SMC round, or per MCMC iteration.
#' @param max_rounds,rho_target SMC: repeat the sweeps (up to `max_rounds`)
#'   until the correlation of particle log-likelihoods before and after falls
#'   below `rho_target`. `max_rounds = 1` gives a fixed amount of work.
#' @param n_sm Split-merge proposals per sweep.
#' @param gamma Dirichlet concentration of the mixture weights.
#' @param p_prior_a,p_prior_b Beta prior on each bloc Yes probability.
#' @param nu_shape,nu_rate Gamma prior on each bloc's dispersion \eqn{\nu_k}
#'   (default mean 100).
#' @param nu_min Hard floor on \eqn{\nu_k}: the Gamma prior is truncated to
#'   \eqn{[\nu_{min}, \infty)}, so no bloc can flatten into a catch-all for
#'   units that fit nowhere else. 0 removes the floor.
#' @param nu_pool,pool_tau,pool_m0,pool_s0 Alternative to the floor: a
#'   hierarchical log-normal prior, \eqn{\log \nu_k \mid \mu \sim N(\mu,
#'   \tau^2)}, \eqn{\mu \sim N(m_0, s_0^2)}. Requires `nu_min = 0`. `NULL`
#'   `pool_m0`/`pool_s0` match the log-scale mean and variance of the Gamma prior.
#' @param p_sd,nu_sd Starting random-walk sds on logit \eqn{p} and log \eqn{\nu}.
#' @param adapt Adapt the random-walk sds (SMC: between steps; MCMC: during burn-in only).
#' @param sm_method Split-merge proposal: `"jain_neal"` (launch-state
#'   restricted Gibbs, the default), `"anchor"` or `"prior"`.
#' @param launch,n_launch,n_p_launch,sm_p_rw_sd,sm_nu_rw_sd,delta_sd,delta_nu_sd
#'   Tuning of the Jain--Neal launch state (see the package source).
#' @param sm_conc0,sm_nu_ref Tuning of the `"anchor"` proposals.
#' @param schedule `"adaptive"` (by conditional ESS) or `"power"`: fixed
#'   \eqn{\beta_r = \zeta (r / n_{rungs})^{a}}.
#' @param n_rungs,gti_pow Used by `schedule = "power"`.
#' @param max_steps Abort the SMC after this many tempering steps.
#' @return A named list of settings.
#' @examples
#' str(bbmm_control())
#' bbmm_control(n_particles = 100, nu_min = 20)$n_particles
#' @export
bbmm_control <- function(n_particles = 200L, cess_target = 0.9, resample_threshold = 0.5, n_sweeps = 2L,
                         max_rounds = 1L, rho_target = 0.5, n_sm = 5L, gamma = 1,
                         p_prior_a = 10, p_prior_b = 10, nu_shape = 4, nu_rate = 0.04, nu_min = 50,
                         nu_pool = FALSE, pool_tau = 0.25, pool_m0 = NULL, pool_s0 = NULL,
                         p_sd = 1, nu_sd = 0.5, adapt = TRUE,
                         sm_method = c("jain_neal", "anchor", "prior"), launch = c("pooled", "mh"),
                         n_launch = 3L, n_p_launch = 1L, sm_p_rw_sd = 0.5, sm_nu_rw_sd = 0.5,
                         delta_sd = 1, delta_nu_sd = 1, sm_conc0 = 20, sm_nu_ref = 100,
                         schedule = c("adaptive", "power"), n_rungs = 100L, gti_pow = 3, max_steps = 100000L) {
  ctl <- as.list(environment())
  ctl$sm_method <- match.arg(sm_method); ctl$launch <- match.arg(launch); ctl$schedule <- match.arg(schedule)
  ctl
}

# Internal: merge user controls into the defaults and build the list the C++ code reads.
cpp_control <- function(K, control, zeta, seed, n_threads, progress, n_particles = NULL) {
  defaults <- bbmm_control()
  if (inherits(control, "list") && length(control)) {
    bad <- setdiff(names(control), names(defaults))
    if (length(bad)) stop("Unknown control setting(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  ctl <- utils::modifyList(defaults, as.list(control))
  if (!is.null(n_particles)) ctl$n_particles <- n_particles
  ctl$K <- as.integer(K)
  if (!(is.numeric(zeta) && length(zeta) == 1 && zeta > 0 && zeta <= 1)) stop("`zeta` must lie in (0, 1].", call. = FALSE)
  if (isTRUE(ctl$nu_pool)) {
    if (ctl$nu_min > 0) stop("nu_pool replaces the floor on nu: set nu_min = 0.", call. = FALSE)
    if (is.null(ctl$pool_m0)) ctl$pool_m0 <- digamma(ctl$nu_shape) - log(ctl$nu_rate)
    if (is.null(ctl$pool_s0)) {
      v <- trigamma(ctl$nu_shape) - ctl$pool_tau^2
      if (v <= 0) stop("pool_tau^2 must be below trigamma(nu_shape) unless pool_s0 is given.", call. = FALSE)
      ctl$pool_s0 <- sqrt(v)
    }
  }
  stopifnot(ctl$K >= 1, ctl$n_particles >= 2, ctl$cess_target > 0, ctl$cess_target < 1,
            ctl$resample_threshold > 0, ctl$resample_threshold <= 1, ctl$n_sweeps >= 1, ctl$n_sm >= 0,
            ctl$gamma > 0, ctl$p_prior_a > 0, ctl$p_prior_b > 0, ctl$nu_shape > 0, ctl$nu_rate > 0,
            ctl$nu_min >= 0)
  ints <- c("n_particles", "n_sweeps", "n_sm", "max_steps", "n_launch", "n_p_launch", "max_rounds")
  ctl[ints] <- lapply(ctl[ints], as.integer)
  cc <- ctl[c("K", ints, "gamma", "p_prior_a", "p_prior_b", "nu_shape", "nu_rate", "cess_target",
              "resample_threshold", "p_sd", "nu_sd", "adapt", "rho_target", "nu_min", "sm_conc0", "sm_nu_ref",
              "sm_p_rw_sd", "sm_nu_rw_sd", "delta_sd", "delta_nu_sd")]
  cc$n_threads <- max(1L, as.integer(n_threads))
  cc$zeta <- zeta
  cc$seed <- as.numeric(seed)
  cc$progress <- isTRUE(progress)
  cc$nu_pool <- isTRUE(ctl$nu_pool)
  if (cc$nu_pool) cc[c("pool_tau", "pool_m0", "pool_s0")] <- ctl[c("pool_tau", "pool_m0", "pool_s0")]
  cc$sm_method <- match(ctl$sm_method, c("prior", "anchor", "jain_neal")) - 1L
  cc$launch_method <- match(ctl$launch, c("mh", "pooled")) - 1L
  cc$beta_grid <- if (ctl$schedule == "power") zeta * ((1:ctl$n_rungs) / ctl$n_rungs)^ctl$gti_pow else numeric(0)
  ctl$zeta <- zeta
  list(ctl = ctl, cc = cc)
}

# Internal: counts as the C++ code wants them (NA for empty cells, integer).
cpp_counts <- function(data) list(X = data$X, M = data$M)

# Internal: an integer seed for the C++ streams, from `seed` or R's stream.
cpp_seed <- function(seed) if (is.null(seed)) sample.int(.Machine$integer.max - 1L, 1L) else as.integer(seed)

#' Marginal likelihood and posterior sample at fixed K by stepping-stone SMC
#'
#' Estimates the evidence \eqn{\log p(x \mid K)} of the Bayesian BBMM with
#' \eqn{K} components by adaptive-tempering sequential Monte Carlo along the
#' path from the prior (\eqn{\beta = 0}) to the tempered posterior
#' (\eqn{\beta = \zeta}): the stepping-stone estimator with resampling. The
#' final weighted particles are a sample from the posterior at that \eqn{K}.
#' Run it for several \eqn{K} with [bbmm_evidence()] to obtain \eqn{p(K \mid x)}.
#'
#' At each step the next \eqn{\beta} is chosen by bisection so that the
#' conditional ESS of the incremental weights is `cess_target`; particles are
#' resampled when the ESS falls below `resample_threshold`, then moved by an
#' MCMC kernel that leaves the current tempered posterior invariant: empty
#' components redrawn from the prior, Jain--Neal split-merge proposals, a
#' collapsed Gibbs sweep of the allocations, random-walk Metropolis updates of
#' logit \eqn{p_{kq}} and log \eqn{\nu_k}. Thermodynamic-integration
#' estimates from the same path are returned for comparison.
#'
#' **Run time.** The cost grows with units x questions and with the number of
#' tempering steps (a few hundred on real data). With the default 200
#' particles, Maine (462 x 72) takes roughly 5--20 minutes per \eqn{K} on 8
#' threads of a recent laptop; a state's precincts (thousands x dozens) can
#' take hours per \eqn{K}. The SMC estimate of \eqn{\log Z} is biased low when
#' particles lag the moving target, so run several replicates if \eqn{K}
#' comparisons are close.
#'
#' @param data A [bbmm_data] object.
#' @param K Number of mixture components.
#' @param zeta Likelihood power in (0, 1].
#' @param control A list of settings; see [bbmm_control()].
#' @param seed Optional integer seed; results are reproducible from it and do
#'   not depend on `n_threads` (on the same platform). If `NULL`, the seed is
#'   drawn from R's random-number stream.
#' @param n_threads Number of threads.
#' @param progress Print a progress line every 25 steps.
#' @param keep_particles Keep the final particles (needed for summaries).
#' @return An object of class `bbmm_ss`: `logZ` (stepping-stone estimate),
#'   `estimates` (SS, TI, corrected TI), `steps` (one row per tempering step),
#'   the final particles `z` (particles x units), `p` (particles x K x Q
#'   array), `nu`, `ll` (untempered complete-data log-likelihood) and
#'   `weights`, `T_dist` (posterior distribution of the number of occupied
#'   components), and the settings.
#' @references Del Moral, P., Doucet, A. and Jasra, A. (2006). Sequential Monte
#'   Carlo samplers. *JRSS B* 68, 411--436.
#'
#'   Xie, W., Lewis, P. O., Fan, Y., Kuo, L. and Chen, M.-H. (2011). Improving
#'   marginal likelihood estimation for Bayesian phylogenetic model selection.
#'   *Systematic Biology* 60, 150--160.
#'
#'   Jain, S. and Neal, R. M. (2004). A split-merge Markov chain Monte Carlo
#'   procedure for the Dirichlet process mixture model. *Journal of
#'   Computational and Graphical Statistics* 13, 158--182.
#' @seealso [bbmm_evidence()], [bbmm_rejuvenate()], [bloc_summary()]
#' @examples
#' tiny <- maine[1:40, 1:6]
#' fit <- bbmm_ss(tiny, K = 2, zeta = 0.3, control = list(n_particles = 50, cess_target = 0.5),
#'                seed = 1, progress = FALSE)
#' fit
#' @export
bbmm_ss <- function(data, K, zeta = 1, control = list(), seed = NULL, n_threads = 1L,
                    progress = interactive(), keep_particles = TRUE) {
  check_data(data)
  sc <- cpp_control(K, control, zeta, cpp_seed(seed), n_threads, progress)
  cm <- cpp_counts(data)
  t0 <- proc.time()[["elapsed"]]
  out <- ss_smc_cpp(cm$X, cm$M, sc$cc)
  out$seconds <- proc.time()[["elapsed"]] - t0
  out$estimates <- c(SS = out$logZ, ti_estimates(out$steps, zeta = zeta,
                     gti_pow = if (sc$ctl$schedule == "power") sc$ctl$gti_pow else numeric(0)))
  out$n_steps <- nrow(out$steps) - 1L
  out$K <- as.integer(K); out$zeta <- zeta; out$control <- sc$ctl; out$seed <- sc$cc$seed
  out$units <- data$units; out$questions <- data$questions
  dimnames(out$p) <- list(NULL, NULL, data$questions)
  T_occ <- apply(out$z, 1, function(z) length(unique(z)))
  out$T_dist <- tapply(out$weights, factor(T_occ, levels = seq_len(K)), sum, default = 0)
  if (!isTRUE(sc$ctl$nu_pool)) out$mu <- NULL
  if (!keep_particles) out[c("z", "p", "nu", "weights", "ll", "mu")] <- NULL
  class(out) <- "bbmm_ss"
  out
}

# Thermodynamic integration from the SMC path: trapezoid (TI), the corrected
# trapezoid of Friel, Hurn & Wyse (2014), and generalised TI on a power grid.
ti_estimates <- function(steps, gti_pow = numeric(0), zeta = 1) {
  b <- steps$beta / zeta; E <- zeta * steps$mean_ll; V <- zeta^2 * steps$var_ll
  db <- diff(b)
  ti <- sum(db * (head(E, -1) + tail(E, -1)) / 2)
  ti_corr <- ti - sum(db^2 / 12 * diff(V))
  gti <- vapply(gti_pow, function(a) {
    t <- b^(1 / a)
    g <- a * t^(a - 1) * E
    g[t == 0] <- if (a > 1) 0 else E[t == 0]
    sum(diff(t) * (head(g, -1) + tail(g, -1)) / 2)
  }, numeric(1))
  if (length(gti_pow)) names(gti) <- paste0("GTI_", gti_pow)
  c(TI = ti, TI_corr = ti_corr, gti)
}

#' @export
print.bbmm_ss <- function(x, ...) {
  cat(sprintf("BBMM stepping-stone SMC: K = %d, zeta = %g, %d particles, %d steps, %.0f s\n",
              x$K, x$zeta, x$control$n_particles, x$n_steps, x$seconds))
  cat(sprintf("log p(x | K) = %.2f (SS); TI %.2f, corrected TI %.2f\n", x$logZ,
              x$estimates[["TI"]], x$estimates[["TI_corr"]]))
  Td <- x$T_dist[x$T_dist > 0.005]
  cat("Posterior on occupied components:", paste(sprintf("T=%s: %.2f", names(Td), Td), collapse = ", "), "\n")
  invisible(x)
}

#' Posterior on the number of blocs: stepping-stone runs over a range of K
#'
#' Runs [bbmm_ss()] for each \eqn{K} (and `n_rep` replicates) and combines
#' the evidences with a prior on \eqn{K} into \eqn{p(K \mid x)}, and with each
#' run's particles into \eqn{p(T \mid x)}, the posterior on the number of
#' *occupied* blocs \eqn{T \le K}.
#'
#' The default prior is geometric with \eqn{r = 0.1} on \eqn{1, \dots, 50},
#' \eqn{P(K) \propto r (1 - r)^{K - 1}}. With replicates, \eqn{p(K \mid x)}
#' uses the mean of the replicates' \eqn{\log Z}; Monte Carlo error in
#' \eqn{\log Z} is not propagated (compare `sd` across replicates with the
#' gaps between \eqn{K}).
#'
#' **Run time.** One [bbmm_ss()] run per \eqn{K} and replicate: see its
#' help. On Maine (462 x 72), \eqn{K = 2, \dots, 8} takes on the order of an
#' hour on 8 threads.
#'
#' @inheritParams bbmm_ss
#' @param K Integer vector of component counts.
#' @param n_rep Replicates per K (different seeds).
#' @param prior_K Function of K returning the log prior; default geometric(0.1) on 1..50.
#' @param keep_fits Keep the fitted runs (with particles) for summaries.
#' @return An object of class `bbmm_evidence`: `runs` (one row per K and
#'   replicate: logZ, TI, steps, seconds), `posterior_K` (K, mean and sd of
#'   logZ, prior, posterior), `posterior_T` (T, posterior) and `fits`.
#' @seealso [plot_evidence()], [bloc_summary()]
#' @examples
#' tiny <- maine[1:40, 1:6]
#' ev <- bbmm_evidence(tiny, K = 1:3, zeta = 0.3, seed = 1, progress = FALSE,
#'                     control = list(n_particles = 50, cess_target = 0.5))
#' ev
#' @export
bbmm_evidence <- function(data, K, zeta = 1, n_rep = 1L, control = list(), prior_K = NULL, seed = NULL,
                          n_threads = 1L, progress = interactive(), keep_fits = TRUE) {
  check_data(data)
  K <- sort(unique(as.integer(K)))
  seeds <- start_seeds(seed, length(K) * n_rep)
  prior_K <- prior_K %||% log_prior_K_geometric
  fits <- list(); rows <- list(); s <- 0L
  for (k in K) for (r in seq_len(n_rep)) {
    s <- s + 1L
    if (progress) message(sprintf("K = %d, replicate %d", k, r))
    f <- bbmm_ss(data, k, zeta, control, seed = seeds[s], n_threads = n_threads, progress = progress,
                 keep_particles = keep_fits)
    rows[[s]] <- data.frame(K = k, rep = r, logZ = f$logZ, TI = f$estimates[["TI"]],
                            TI_corr = f$estimates[["TI_corr"]], n_steps = f$n_steps, seconds = f$seconds)
    if (keep_fits) fits[[s]] <- f
  }
  runs <- do.call(rbind, rows)
  out <- list(runs = runs, fits = if (keep_fits) fits, zeta = zeta, prior_K = prior_K)
  class(out) <- "bbmm_evidence"
  evidence_posteriors(out)
}

# Internal: recompute posterior_K and posterior_T from runs (and fits).
evidence_posteriors <- function(ev) {
  runs <- ev$runs
  agg <- do.call(rbind, lapply(split(runs, runs$K), function(t)
    data.frame(K = t$K[1], logZ = mean(t$logZ), sd = if (nrow(t) > 1) sd(t$logZ) else NA_real_, n_rep = nrow(t))))
  agg$log_prior <- ev$prior_K(agg$K)
  agg$posterior <- posterior_K(agg$logZ, agg$K, agg$log_prior)
  ev$posterior_K <- agg
  if (!is.null(ev$fits)) {
    Kmax <- max(runs$K)
    pT <- numeric(Kmax)
    for (s in seq_len(nrow(runs))) {
      f <- ev$fits[[s]]
      w <- agg$posterior[agg$K == f$K] / agg$n_rep[agg$K == f$K]
      pT[seq_along(f$T_dist)] <- pT[seq_along(f$T_dist)] + w * as.numeric(f$T_dist)
    }
    ev$posterior_T <- data.frame(T = seq_len(Kmax), posterior = pT)
  }
  ev
}

#' Posterior over K from log evidences
#'
#' @param logZ Log marginal likelihoods, one per K.
#' @param K The component counts.
#' @param log_prior Log prior probabilities of K (default: uniform over the K given).
#' @return Named vector of posterior probabilities.
#' @examples
#' posterior_K(c(-1002, -1000, -1001), K = 2:4, log_prior = log_prior_K_geometric(2:4))
#' @export
posterior_K <- function(logZ, K, log_prior = rep(0, length(K))) {
  lp <- logZ + log_prior
  w <- exp(lp - max(lp))
  stats::setNames(w / sum(w), K)
}

#' Geometric prior on the number of components
#'
#' \eqn{P(K) = r (1 - r)^{K - 1}} truncated to \eqn{1, \dots, K_{max}} and renormalised.
#'
#' @param K Component counts.
#' @param r Success probability.
#' @param K_max Truncation point.
#' @return Log prior probabilities.
#' @examples
#' exp(log_prior_K_geometric(1:5))
#' @export
log_prior_K_geometric <- function(K, r = 0.1, K_max = 50) {
  ifelse(K >= 1 & K <= K_max, log(r) + (K - 1) * log1p(-r) - log1p(-(1 - r)^K_max), -Inf)
}

#' @export
print.bbmm_evidence <- function(x, ...) {
  cat(sprintf("BBMM stepping-stone evidence, zeta = %g: %d run(s)\n", x$zeta, nrow(x$runs)))
  pk <- x$posterior_K
  print(data.frame(K = pk$K, logZ = round(pk$logZ, 2), sd = round(pk$sd, 2), posterior = signif(pk$posterior, 3)),
        row.names = FALSE)
  if (!is.null(x$posterior_T)) {
    pt <- x$posterior_T[x$posterior_T$posterior > 0.005, ]
    cat("Posterior on occupied blocs T:", paste(sprintf("%d: %.2f", pt$T, pt$posterior), collapse = ", "), "\n")
  }
  invisible(x)
}

#' Enlarge and refresh a stepping-stone run's particle sample
#'
#' Resamples `n_final` particles from a finished [bbmm_ss()] run and moves
#' each one `n_rounds` times with the run's final kernel, so that duplicate
#' particles separate. Useful before summaries: SMC ends with many copies of
#' a few particles. The evidence is unchanged.
#'
#' @param fit A `bbmm_ss` object with particles.
#' @param data The [bbmm_data] the run used.
#' @param n_final Number of particles after rejuvenation.
#' @param n_rounds Kernel rounds per particle.
#' @param n_threads Number of threads.
#' @return The `bbmm_ss` object with new particles (equal weights) and a
#'   `rejuvenation` entry.
#' @examples
#' tiny <- maine[1:40, 1:6]
#' fit <- bbmm_ss(tiny, K = 2, zeta = 0.3, control = list(n_particles = 50, cess_target = 0.5),
#'                seed = 1, progress = FALSE)
#' fit2 <- bbmm_rejuvenate(fit, tiny, n_final = 100, n_rounds = 2)
#' @export
bbmm_rejuvenate <- function(fit, data, n_final = 1000L, n_rounds = 10L, n_threads = 1L) {
  stopifnot(inherits(fit, "bbmm_ss"), !is.null(fit$z))
  check_data(data)
  ctl <- fit$control
  sc <- cpp_control(fit$K, ctl[intersect(names(ctl), names(bbmm_control()))], fit$zeta, fit$seed, n_threads, FALSE)
  cm <- cpp_counts(data)
  last <- fit$steps[nrow(fit$steps), ]
  t0 <- proc.time()[["elapsed"]]
  r <- ss_rejuvenate_cpp(cm$X, cm$M, sc$cc, fit$zeta, fit$z, fit$p, fit$nu, fit$weights,
                         as.integer(n_final), as.integer(n_rounds), last$p_sd, last$nu_sd,
                         if (is.null(fit$mu)) numeric(0) else fit$mu)
  fit[c("z", "p", "nu", "ll", "weights")] <- r[c("z", "p", "nu", "ll", "weights")]
  dimnames(fit$p) <- list(NULL, NULL, fit$questions)
  if (isTRUE(ctl$nu_pool)) fit$mu <- r$mu
  T_occ <- apply(fit$z, 1, function(z) length(unique(z)))
  fit$T_dist <- tapply(fit$weights, factor(T_occ, levels = seq_len(fit$K)), sum, default = 0)
  fit$rejuvenation <- list(n_final = n_final, n_rounds = n_rounds, n_distinct_sources = length(unique(r$source)),
                           acceptance = r$acceptance, seconds = proc.time()[["elapsed"]] - t0)
  fit
}

#' Combine separately run stepping-stone fits into a posterior over K
#'
#' Builds the same object as [bbmm_evidence()] from a list of [bbmm_ss()]
#' fits, e.g. one per K run as separate jobs on a cluster. Fits with the same
#' K count as replicates.
#'
#' @param fits A list of `bbmm_ss` objects on the same data and zeta.
#' @param prior_K Function of K returning the log prior; default [log_prior_K_geometric()].
#' @return A `bbmm_evidence` object.
#' @examples
#' tiny <- maine[1:40, 1:6]
#' ctl <- list(n_particles = 50, cess_target = 0.5)
#' fits <- lapply(1:2, function(K)
#'   bbmm_ss(tiny, K, zeta = 0.3, control = ctl, seed = K, progress = FALSE))
#' evidence_from_fits(fits)
#' @export
evidence_from_fits <- function(fits, prior_K = NULL) {
  if (!all(vapply(fits, inherits, TRUE, "bbmm_ss"))) stop("`fits` must be a list of bbmm_ss objects.", call. = FALSE)
  zeta <- unique(vapply(fits, `[[`, 1, "zeta"))
  if (length(zeta) != 1) stop("The fits use different values of zeta.", call. = FALSE)
  if (length(unique(lapply(fits, `[[`, "units"))) != 1) stop("The fits use different data.", call. = FALSE)
  fits <- fits[order(vapply(fits, `[[`, 1L, "K"))]
  Ks <- vapply(fits, `[[`, 1L, "K")
  runs <- data.frame(K = Ks, rep = stats::ave(Ks, Ks, FUN = seq_along), logZ = vapply(fits, `[[`, 1, "logZ"),
                     TI = vapply(fits, function(f) f$estimates[["TI"]], 1),
                     TI_corr = vapply(fits, function(f) f$estimates[["TI_corr"]], 1),
                     n_steps = vapply(fits, `[[`, 1L, "n_steps"), seconds = vapply(fits, `[[`, 1, "seconds"))
  keep <- !vapply(fits, function(f) is.null(f$z), TRUE)
  out <- structure(list(runs = runs, fits = if (all(keep)) fits, zeta = zeta,
                        prior_K = prior_K %||% log_prior_K_geometric), class = "bbmm_evidence")
  evidence_posteriors(out)
}
