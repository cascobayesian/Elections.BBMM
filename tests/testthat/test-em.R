test_that("EM log-likelihood matches a direct beta-binomial evaluation", {
  d <- maine[1:30, 1:5]
  x <- d$X; n <- d$M; x[is.na(x)] <- 0; n[is.na(n)] <- 0; storage.mode(x) <- storage.mode(n) <- "double"
  p <- matrix(c(0.3, 0.5, 0.7, 0.4, 0.6, 0.35, 0.55, 0.45, 0.62, 0.5), 2, 5)
  nu <- 80
  direct <- sapply(1:2, function(k) rowSums(sapply(1:5, function(q) {
    a <- nu * p[k, q]; b <- nu * (1 - p[k, q])
    lchoose(n[, q], x[, q]) + lbeta(x[, q] + a, n[, q] - x[, q] + b) - lbeta(a, b)
  })))
  expect_equal(Elections.BBMM:::bb_loglik_em_cpp(x, n, p, nu), direct, tolerance = 1e-9, ignore_attr = TRUE)
})

test_that("EM is reproducible, does not disturb the RNG, and selects K", {
  d <- maine[1:80, 1:12]
  set.seed(5); before <- runif(1); set.seed(5)
  f1 <- bbmm_em(d, K = 2, zeta = 0.3, n_starts = 2, seed = 11)
  expect_identical(runif(1), before)
  f2 <- bbmm_em(d, K = 2, zeta = 0.3, n_starts = 2, seed = 11)
  expect_identical(f1$loglik, f2$loglik)
  expect_true(f1$converged)
  expect_equal(rowSums(f1$gamma), rep(1, 80), ignore_attr = TRUE)
  sel <- bbmm_em_select(d, K = 1:3, zeta = 0.3, n_starts = 2, seed = 1)
  expect_identical(sel$table$K, 1:3)
  expect_true(all(diff(sel$table$loglik) > -1e-6))
  expect_named(sel$selected, c("BIC", "min_share_3%", "min_share_5%"))
})
