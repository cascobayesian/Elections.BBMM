test_that("bbmm_data validates and stores counts", {
  d <- bbmm_data(matrix(c(3, 0, 5, 2), 2), matrix(c(10, 0, 9, 4), 2), units = c("a", "b"))
  expect_s3_class(d, "bbmm_data")
  expect_true(is.na(d$X[2, 1]) && is.na(d$M[2, 1]))
  expect_identical(d$questions, c("Q1", "Q2"))
  expect_error(bbmm_data(matrix(c(11, 1), 2), matrix(c(10, 2), 2)), "0 <= yes <= total")
  expect_error(bbmm_data(matrix(1:2, 2), matrix(1:3, 3)), "same dimensions")
})

test_that("wide tables and subsetting work", {
  df <- data.frame(id = c("A", "B", "C"), Q1_Y = c(3, 4, 5), Q1_N = c(1, 2, 3), Q2_Y = c(0, 1, 2), Q2_N = c(5, 5, 5))
  d <- bbmm_data_from_wide(df, "id")
  expect_identical(unname(d$M[, "Q2"]), c(5L, 6L, 7L))
  s <- d[c("C", "A"), "Q1"]
  expect_identical(s$units, c("C", "A"))
  expect_identical(dim(s), c(2L, 1L))
})

test_that("the Maine data are complete", {
  expect_identical(dim(maine), c(462L, 72L))
  expect_identical(maine$questions, maine_questions$question_id)
  expect_true(all(maine$X <= maine$M, na.rm = TRUE))
})
