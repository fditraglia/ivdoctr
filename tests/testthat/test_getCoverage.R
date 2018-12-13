context("Interval Coverage")

test_that("getCoverage properly computes coverage", {
  data <- data.frame(min = c(1, 2, 3, 2), max = c(5, 7, 8, 6))
  guess <- c(2, 7)
  expected_answer <- 0.5
  ans <- getCoverage(data, guess)
  expect_equal(expected_answer, ans)
})

test_that("getCoverage properly computes coverage", {
  data <- data.frame(min = c(1, 2, 3, 2), max = c(5, 7, 8, 6))
  guess <- c(0, 0)
  expected_answer <- 0
  ans <- getCoverage(data, guess)
  expect_equal(expected_answer, ans)
})

test_that("getInterval properly computes covering interval", {
  data <- data.frame(min = c(1, 2, 3, 2), max = c(5, 7, 8, 6))
  expected_answer <- c(2, 7)
  ans <- getInterval(data, center = c(4.5, 4.5), conf = 0.5)
  expect_equal(expected_answer, ans, tol = 1e-6)
})

test_that("getInterval properly computes covering interval", {
  data <- data.frame(min = c(1, 2, 3, 2), max = c(5, 7, 8, 6))
  expected_answer <- c(1, 8)
  ans <- getInterval(data, c(4.5, 4.5), conf = 0.95)
  expect_equal(expected_answer, ans, tol = 1e-6)
})

