context("Appendix Functions")

test_that("g_functionA2 properly computes value", {
  kappa <- 1
  r_TstarU <- 0
  obs_draws <- data.frame(r_Tz = 1,
                          r_Ty = 1,
                          r_zy = 1,
                          s2_y = 1,
                          s2_z = 1,
                          s_Tz = 1)
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
                          s2_y = 0.5,
                          s2_z = 1,
                          s_Tz = 1)
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
                          s2_y = 1,
                          s2_z = 1,
                          s_Tz = 1)
  expected_answer <- -1
  ans <- g_functionA2(kappa, r_TstarU, obs_draws)
  expect_equal(expected_answer, ans)
})

test_that("g_functionA2 properly computes value (Vectorized over kappas)", {
  kappa <- rep(1, 10)
  r_TstarU <- 0
  obs_draws <- data.frame(r_Tz = 1,
                          r_Ty = 1,
                          r_zy = 1,
                          s2_y = 1,
                          s2_z = 1,
                          s_Tz = 1)
  expected_answer <- rep(0, 10)
  ans <- g_functionA2(kappa, r_TstarU, obs_draws)
  expect_equal(expected_answer, ans)
})

test_that("b_functionA3 properly computes value", {
  g <- 0
  psi <- -1
  obs_draws <- data.frame(s_Tz = 1,
                          s_zy = 1)
  expected_answer <- cbind(beta_lower = 0, beta_upper = 0)
  ans <- b_functionA3(obs_draws, g, psi)
  expect_equal(expected_answer, ans)
})

test_that("b_functionA3 properly computes value (Vectorized)", {
  g <- rep(0, 10)
  psi <- rep(-1, 10)
  obs_draws <- data.frame(s_Tz = 1,
                          s_zy = 1)
  expected_answer <- cbind(beta_lower = 0, beta_upper = 0)
  ans <- b_functionA3(obs_draws, g, psi)
  expect_equal(expected_answer, ans)
})

test_that("b_functionA3 properly computes value (Vectorized)", {
  g <- rep(0, 11)
  psi <- seq(0, 1, by = 0.1)
  obs_draws <- data.frame(s_Tz = 1,
                          s_zy = 1)
  expected_answer <- cbind(beta_lower = 1, beta_upper = 2)
  ans <- b_functionA3(obs_draws, g, psi)
  expect_equal(expected_answer, ans)
})
