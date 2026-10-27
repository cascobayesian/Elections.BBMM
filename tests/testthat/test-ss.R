# Exact evidence by enumeration on a tiny data set: K = 1, so the evidence is a
# product over questions of a one-dimensional integral over p, and nu is
# integrated by quadrature over its truncated Gamma prior.
test_that("stepping-stone evidence matches the exact value for K = 1", {
  X <- matrix(c(6, 3, 8, 5, 2, 7), 3); M <- matrix(c(10, 9, 12, 8, 6, 11), 3)
  d <- bbmm_data(X, M)
  zeta <- 0.5; nu_min <- 50
  lik_given_nu <- function(nu) prod(sapply(1:2, function(q) stats::integrate(function(p) sapply(p, function(pp) {
    a <- nu * pp; b <- nu * (1 - pp)
    exp(zeta * sum(lbeta(X[, q] + a, M[, q] - X[, q] + b) - lbeta(a, b))) * stats::dbeta(pp, 10, 10)
  }), 0, 1)$value))
  dens <- function(nu) dgamma(nu, 4, rate = 0.04) / pgamma(nu_min, 4, rate = 0.04, lower.tail = FALSE)
  exact <- log(stats::integrate(function(v) sapply(v, function(nu) lik_given_nu(nu) * dens(nu)), nu_min, Inf)$value)
  est <- replicate(3, bbmm_ss(d, K = 1, zeta = zeta, seed = sample.int(1e6, 1), progress = FALSE,
                              control = list(n_particles = 400))$logZ)
  expect_lt(abs(mean(est) - exact), 0.1)
})

test_that("SS results do not depend on the number of threads", {
  d <- maine[1:40, 1:6]
  ctl <- list(n_particles = 40, cess_target = 0.5)
  a <- bbmm_ss(d, K = 2, zeta = 0.3, control = ctl, seed = 3, n_threads = 1, progress = FALSE)
  b <- bbmm_ss(d, K = 2, zeta = 0.3, control = ctl, seed = 3, n_threads = 2, progress = FALSE)
  expect_identical(a$logZ, b$logZ)
  expect_identical(a$z, b$z)
  expect_equal(sum(a$T_dist), 1)
  ev <- bbmm_evidence(d, K = 1:2, zeta = 0.3, control = ctl, seed = 1, progress = FALSE)
  expect_equal(sum(ev$posterior_K$posterior), 1)
  expect_equal(sum(ev$posterior_T$posterior), 1)
  r <- bbmm_rejuvenate(a, d, n_final = 60, n_rounds = 1)
  expect_identical(nrow(r$z), 60L)
})

test_that("posterior_K and the geometric prior behave", {
  expect_equal(sum(exp(log_prior_K_geometric(1:50))), 1)
  expect_equal(unname(posterior_K(c(0, log(3)), 1:2)), c(0.25, 0.75))
})
