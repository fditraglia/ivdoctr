
test_that("candidate1 computes properly 1", {
  obs <- list(r_Ty = 0,
              r_Tz = 0,
              r_zy = 0)
  r_TstarU_upper <- 1
  r_TstarU_lower <- 0
  k_upper <- 1
  k_lower <- 0.1
  expected_answer <- c(0, 0)
  ans <- candidate1(r_TstarU_upper, r_TstarU_lower, k_upper, k_lower, obs)
  expect_equal(expected_answer, ans)
})

test_that("candidate1 is appropriately vectorized wrt r_TstarU and k", {
  obs <- list(r_Ty = 0,
              r_Tz = 0,
              r_zy = 0)
  r_TstarU_upper <- c(1, 1)
  r_TstarU_lower <- c(0, 0)
  k_upper <- c(1, 1)
  k_lower <- c(0.1, 0.1)
  expected_answer <- cbind(min_corner = c(0, 0), max_corner = c(0, 0))
  ans <- candidate1(r_TstarU_upper, r_TstarU_lower, k_upper, k_lower, obs)
  expect_equal(expected_answer, ans)
})

test_that("candidate2 properly computes solutions", {
  obs <- list(r_Ty = 0,
              r_Tz = 0,
              r_zy = 1)
  r_TstarU_upper <- 1
  r_TstarU_lower <- 0
  k_upper <- 1
  k_lower <- 0.1
  expected_answer <- cbind(1, 1)
  ans <- candidate2(r_TstarU_upper, r_TstarU_lower, k_upper, k_lower, obs)
  expect_equivalent(expected_answer, ans)
})

test_that("candidate2 is vectorized wrt r_TstarU and k", {
  obs <- list(r_Ty = 0,
              r_Tz = 0,
              r_zy = 1)
  r_TstarU_upper <- c(1, 1)
  r_TstarU_lower <- c(0, 0)
  k_upper <- c(1, 1)
  k_lower <- c(0.1, 0.1)
  expected_answer <- cbind(min_edge = c(1, 1), max_edge = c(1, 1))
  ans <- candidate2(r_TstarU_upper, r_TstarU_lower, k_upper, k_lower, obs)
  expect_equivalent(expected_answer, ans)
})
