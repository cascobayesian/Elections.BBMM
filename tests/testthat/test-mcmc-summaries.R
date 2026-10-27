test_that("fixed-K MCMC runs, is reproducible, and summarises", {
  d <- maine[1:50, 1:8]
  a <- bbmm_mcmc(d, K = 3, zeta = 0.3, n_iter = 40, n_burn = 20, thin = 2, n_chains = 2, seed = 4, progress = FALSE)
  b <- bbmm_mcmc(d, K = 3, zeta = 0.3, n_iter = 40, n_burn = 20, thin = 2, n_chains = 2, seed = 4, progress = FALSE,
                 n_threads = 2)
  expect_identical(a$z, b$z)
  expect_identical(dim(a$z), c(40L, 50L))
  expect_identical(dim(a$p), c(40L, 3L, 8L))
  expect_true(all(a$nu >= 50))
  s <- bloc_summary(a, d)
  expect_s3_class(s, "bloc_summary")
  expect_equal(diag(s$psm), rep(1, 50), ignore_attr = TRUE)
  expect_true(all(s$assignment$certainty >= 0 & s$assignment$certainty <= 1))
  expect_identical(sum(s$blocs$n_units), 50L)
  expect_true(all(diff(s$blocs$mean_yes) <= 0))
  pr <- bbmm_mcmc(d, K = 2, zeta = 0.3, n_iter = 10, n_burn = 5, init = "prior", seed = 1, progress = FALSE)
  expect_identical(nrow(pr$z), 20L)
})

test_that("EM summaries, distinctive questions and question fit", {
  d <- maine[1:80, 1:12]
  s <- bloc_summary(bbmm_em(d, K = 3, zeta = 0.3, n_starts = 2, seed = 2), d)
  expect_null(s$psm)
  expect_equal(rowSums(s$membership), rep(1, 80), tolerance = 1e-6, ignore_attr = TRUE)
  dq <- distinctive_questions(s, 2)
  expect_true(all(table(dq$bloc) <= 2))
  qf <- question_fit(s, d)
  expect_true(all(qf$r2 <= 1))
  expect_error(bloc_summary(bbmm_em(d, K = 2, n_starts = 1, seed = 1), maine[1:80, 1:11]), "does not match")
})

test_that("partition helpers", {
  expect_equal(adjusted_rand(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)
  expect_identical(Elections.BBMM:::match_labels(c(3, 3, 1, 1), c(1, 1, 2, 2), 2), c(3L, 1L))
})

test_that("plots build", {
  skip_if_not_installed("ggplot2")
  d <- maine[1:60, 1:10]
  em <- bbmm_em(d, K = 3, zeta = 0.3, n_starts = 2, seed = 1)
  s <- bloc_summary(em, d)
  expect_s3_class(plot_bloc_profiles(s), "ggplot")
  expect_s3_class(plot_js_pca(d, s), "ggplot")
  mc <- bbmm_mcmc(d, K = 3, zeta = 0.3, n_iter = 20, n_burn = 10, seed = 1, progress = FALSE, init = em)
  expect_s3_class(plot_trace(mc), "ggplot")
  expect_s3_class(plot_coclustering(bloc_summary(mc, d)), "ggplot")
  skip_if_not_installed("sf")
  expect_s3_class(plot_bloc_map(s, maine_towns()), "ggplot")
})
