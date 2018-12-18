context("Psi Bounds")

test_that("get_psi_upper properly computes upper bound", {
  p <- 0.5
  kappa <- 0.5
  s2_T <- 1
  expected_answer <- -1
  ans <- get_psi_upper(s2_T, p, kappa)
  expect_equal(expected_answer, ans)
})

test_that("get_psi_upper properly computes upper bound (Vectorized)", {
  p <- 0.5
  kappa <- rep(0.5, 10)
  s2_T <- rep(1, 10)
  expected_answer <- rep(-1, 10)
  ans <- get_psi_upper(s2_T, p, kappa)
  expect_equal(expected_answer, ans)
})

test_that("get_psi_upper properly computes upper bound (Vectorized)", {
  p <- 0.25
  kappa <- rep(0.5, 10)
  s2_T <- rep(1, 10)
  expected_answer <- rep(-2 / 3, 10)
  ans <- get_psi_upper(s2_T, p, kappa)
  expect_equal(expected_answer, ans)
})

test_that("get_psi_lower properly computes lower bound", {
  p <- 0.5
  kappa <- 1
  s2_T <- 1
  expected_answer <- 0
  ans <- get_psi_lower(s2_T, p, kappa)
  expect_equal(expected_answer, ans)
})

test_that("get_psi_lower properly computes lower bound (Vectorized)", {
  p <- 0.5
  kappa <- rep(1, 10)
  s2_T <- rep(1, 10)
  expected_answer <- rep(0, 10)
  ans <- get_psi_lower(s2_T, p, kappa)
  expect_equal(expected_answer, ans)
})

test_that("get_psi_lower properly computes lower bound (Vectorized)", {
  p <- 0.5
  kappa <- rep(0.91, 10)
  s2_T <- rep(1, 10)
  expected_answer <- rep(-0.2, 10)
  ans <- get_psi_lower(s2_T, p, kappa)
  expect_equal(expected_answer, ans)
})

test_that("get_L properly computes L", {
  draws <- data.frame(r_Ty = 1,
                      r_Tz = 1,
                      r_zy = 0)
  expected_answer <- 2
  ans <- get_L(draws)
  expect_equal(expected_answer, ans)
})

test_that("get_L properly computes L", {
  draws <- data.frame(r_Ty = rep(1, 10),
                      r_Tz = rep(1, 10),
                      r_zy = rep(0, 10))
  expected_answer <- rep(2, 10)
  ans <- get_L(draws)
  expect_equal(expected_answer, ans)
})

test_that("get_alpha_bounds properly computes bounds", {
  p <- 0.5
  draws <- data.frame(s2_T = 1,
                      r_Ty = 1,
                      r_Tz = 1,
                      r_zy = 0)
  expected_answer <- list(a0 = -2, a1 = -2)
  ans <- get_alpha_bounds(draws, p)
  expect_equal(expected_answer, ans)
})
