context("Appendix Functions")

test_that("g_functionA2 properly computes value", {
  kappa <- 1
  r_TstarU <- 0
  obs_draws <- data.frame(r_Tz = 1,
                          r_Ty = 1,
                          r_zy = 1,
                          s2_y = 1)
  expected_answer <- 0
  ans <- g_functionA2(kappa, r_TstarU, obs_draws)
  expect_equal(expected_answer, ans)
})

test_that("g_functionA2 properly computes value", {
  kappa <- 1
  r_TstarU <- 0
  obs_draws <- data.frame(r_Tz = 1,
                          r_Ty = 1,
                          r_zy = 1,
                          s2_y = 0.5)
  expected_answer <- 0
  ans <- g_functionA2(kappa, r_TstarU, obs_draws)
  expect_equal(expected_answer, ans)
})

test_that("g_functionA2 properly computes value", {
  kappa <- 1
  r_TstarU <- 0
  obs_draws <- data.frame(r_Tz = 1,
                          r_Ty = 1,
                          r_zy = 0,
                          s2_y = 1)
  expected_answer <- -1
  ans <- g_functionA2(kappa, r_TstarU, obs_draws)
  expect_equal(expected_answer, ans)
})

test_that("g_functionA2 properly computes value (Vectorized over obs_draws)", {
  kappa <- 1
  r_TstarU <- 0
  obs_draws <- data.frame(r_Tz = 1,
                          r_Ty = 1,
                          r_zy = 1,
                          s2_y = 1)
  expected_answer <- 0
  ans <- g_functionA2(kappa, r_TstarU, obs_draws)
  expect_equal(expected_answer, ans)
})
