context("Identified Set")

test_that("Projecting out without controls returns same data", {
  y <- c(2, 4, 6, 8)
  x <- 2 * y
  z <- seq.int(1:4)
  Tobs <- c(3, 7, 5, 2)
  Data <- data.frame(y, x, z, Tobs)
  Sigma <- cov(Data)
  Rho <- cov2cor(Sigma)
  expected_answer <- list(n = 4,
                          T_Rsq = 0,
                          z_Rsq = 0,
                          s2_T = Sigma["Tobs", "Tobs"],
                          s2_y = Sigma["y", "y"],
                          s2_z = Sigma["z", "z"],
                          s_Ty = Sigma["Tobs", "y"],
                          s_Tz = Sigma["Tobs", "z"],
                          s_zy = Sigma["z", "y"],
                          r_Ty = Rho["Tobs", "y"],
                          r_Tz = Rho["Tobs", "z"],
                          r_zy = Rho["z", "y"])
  ans <- get_observables(y_name = "y", T_name = "Tobs", z_name = "z",
                         data = Data, controls = NULL)
  expect_equal(expected_answer, ans)
})

test_that("Projecting out with controls returns projected out data", {
  y <- c(2, 4, 6, 8)
  x <- 2 * y
  z <- seq.int(1:4)
  Tobs <- c(3, 7, 5, 2)
  data <- data.frame(y, x, z, Tobs)
  newY <- resid(lm(data = data, y ~ x))
  newZ <- resid(lm(data = data, z ~ x))
  newT <- resid(lm(data = data, Tobs ~ x))
  newData <- data.frame(newY, newZ, newT)
  sigma <- cov(newData)
  rho <- suppressWarnings(cov2cor(sigma)) # Data is constructed to be a perfect fit.
  reg <- lm(data = data, Tobs ~ x)
  expected_answer <- list(n = 4,
                          T_Rsq = summary(reg)$r.squared,
                          z_Rsq = 1,
                          s2_T = sigma["newT", "newT"],
                          s2_y = sigma["newY", "newY"],
                          s2_z = sigma["newZ", "newZ"],
                          s_Ty = sigma["newT", "newY"],
                          s_Tz = sigma["newT", "newZ"],
                          s_zy = sigma["newZ", "newY"],
                          r_Ty = rho["newT", "newY"],
                          r_Tz = rho["newT", "newZ"],
                          r_zy = rho["newZ", "newY"])
  ans <- suppressWarnings(get_observables(y_name = "y", T_name = "Tobs",
                                          z_name = "z", data = data, controls = "x"))
  expect_equal(expected_answer, ans)
})

test_that("get_k_bounds_unrest properly computes bounds 1", {
  obs <- list(r_Ty = 0.5,
              r_Tz = 0.5,
              r_zy = 0.5,
              T_Rsq = 0)
  ans <- get_k_bounds_unrest(obs, tilde = TRUE)
  expected_answer <- list(Lower = 1 / 3, Upper = 1)
  expect_equal(expected_answer, ans)
})

test_that("get_k_bounds_unrest properly computes bounds 2", {
  obs <- list(r_Ty = 0.5,
              r_Tz = 0.5,
              r_zy = 0.5,
              T_Rsq = 0)
  ans <- get_k_bounds_unrest(obs, tilde = FALSE)
  expected_answer <- list(Lower = 1 / 3, Upper = 1)
  expect_equal(expected_answer, ans)
})

test_that("get_k_bounds_unrest properly computes bounds 3", {
  obs <- list(r_Ty = 0.5,
              r_Tz = 0.5,
              r_zy = 0.5,
              T_Rsq = 0.5)
  ans <- get_k_bounds_unrest(obs, tilde = TRUE)
  expected_answer <- list(Lower = 1 / 3, Upper = 1)
  expect_equal(expected_answer, ans)
})

test_that("get_k_bounds_unrest properly computes bounds 4", {
  obs <- list(r_Ty = 0.5,
              r_Tz = 0.5,
              r_zy = 0.5,
              T_Rsq = 0.5)
  ans <- get_k_bounds_unrest(obs, tilde = FALSE)
  expected_answer <- list(Lower = 2 / 3, Upper = 1)
  expect_equal(expected_answer, ans)
})

test_that("get_k_bounds_unrest properly computes bounds 5. This tests for
          vectorization as well", {
    set.seed(1234)
    nsim <- 500
    lengthGrid <- 500

    # Generating a bunch of correlation draws
    sims_rho_Tz <- 2 * runif(nsim) - 1
    sims_rho_Ty <- 2 * runif(nsim) - 1
    sims_rho_zy <- 2 * runif(nsim) - 1
    obs <- list(r_Ty = sims_rho_Ty,
                r_zy = sims_rho_zy,
                r_Tz = sims_rho_Tz,
                T_Rsq = rep(0, nsim))

    # Getting unrestricted bounds for kappa
    k_bounds <- get_k_bounds_unrest(obs, tilde = FALSE)
    k_lower <- ((sims_rho_Ty^2) + (sims_rho_Tz^2) - 2 * sims_rho_Ty * sims_rho_Tz *
                  sims_rho_zy) / (1 - (sims_rho_zy^2))
    k_upper <- rep(1, nsim)

    expect_equal(k_bounds$Lower, k_lower)
    expect_equal(k_bounds$Upper, k_upper)
})

test_that("get_r_uz_bounds_unrest properly computes bounds 1", {
  obs <- list(r_Ty = 0.5,
              r_Tz = 0.5,
              r_zy = 0.5,
              T_Rsq = 0.5)
  ans <- get_r_uz_bounds_unrest(obs)
  expected_answer <- list(Lower = -1, Upper = 0.5 / sqrt(1 / 3))
  expect_equal(expected_answer, ans)
})

test_that("get_r_uz_bounds_unrest properly computes bounds 2", {
  obs <- list(r_Ty = -0.5,
              r_Tz = 0.5,
              r_zy = 0.5,
              T_Rsq = 0.5)
  ans <- get_r_uz_bounds_unrest(obs)
  expected_answer <- list(Lower = -0.5, Upper = 1)
  expect_equal(expected_answer, ans)
})

test_that("get_r_uz_bounds_unrest properly is vectorized", {
  obs <- list(r_Ty = rep(-0.5, 10),
              r_Tz = rep(0.5, 10),
              r_zy = rep(0.5, 10),
              T_Rsq = rep(0.5, 10))
  ans <- get_r_uz_bounds_unrest(obs)
  expected_answer <- list(Lower = rep(-0.5, 10), Upper = rep(1, 10))
  expect_equal(expected_answer, ans)
})

test_that("get_r_TstarU_bounds_unrest properly is vectorized", {
  obs <- list(r_Ty = rep(-0.5, 10),
              r_Tz = rep(0.5, 10),
              r_zy = rep(0.5, 10),
              T_Rsq = rep(0.5, 10))
  ans <- get_r_TstarU_bounds_unrest(obs)
  expected_answer <- list(Lower = rep(-1, 10), Upper = rep(1, 10))
  expect_equal(expected_answer, ans)
})

test_that("get_r_uz computes properly 1", {
  obs <- list(r_Ty = -0.5,
              r_Tz = 0.5,
              r_zy = 0.5)
  r_TstarU <- 1
  k <- 1
  expected_answer <- 0.5
  ans <- get_r_uz(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_r_uz computes properly 2", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 0,
              r_zy = 0)
  r_TstarU <- 0
  k <- 1
  expected_answer <- 0
  ans <- get_r_uz(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_r_uz computes properly 3", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 0,
              r_zy = 1)
  r_TstarU <- 0
  k <- 1
  expected_answer <- 2
  ans <- get_r_uz(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_r_uz handles multiple sigma draws for the same k and r_TstarU", {
  obs <- list(r_Ty = c(sqrt(3/4), sqrt(3/4)),
              r_Tz = c(0, 0),
              r_zy = c(1, 1))
  r_TstarU <- 0
  k <- 1
  expected_answer <- c(2, 2)
  ans <- get_r_uz(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_r_uz handles multiple sigma draws each with a corresponding
          k and r_TstarU", {
  obs <- list(r_Ty = c(sqrt(3/4), sqrt(3/4)),
              r_Tz = c(0, 0),
              r_zy = c(1, 1))
  r_TstarU <- c(0, 0)
  k <- c(1, 1)
  expected_answer <- c(2, 2)
  ans <- get_r_uz(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_r_uz handles multiple r_TstarU and k for the same sigma draws", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 0,
              r_zy = 1)
  r_TstarU <- c(0, 0)
  k <- c(1, 1)
  expected_answer <- c(2, 2)
  ans <- get_r_uz(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_r_uz throws error with values of k and r_TstarU that are the
          wrong dimension", {
  obs <- list(r_Ty = c(sqrt(3/4), sqrt(3/4), sqrt(3/4)),
              r_Tz = c(0, 0, 0),
              r_zy = c(1, 1, 1))
  r_TstarU <- c(0, 0)
  k <- c(1, 1)
  expect_error(get_r_uz(r_TstarU, k, obs))
})

test_that("get_M properly computes M", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 0,
              r_zy = 1)
  r_TstarU <- 0
  k <- 1
  expected_answer <- sqrt(10)
  ans <- get_M(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_M properly computes M for multiple sims but one value of r_TstarU
          and k", {
  obs <- list(r_Ty = c(sqrt(3/4), sqrt(3/4)),
              r_Tz = c(0, 0),
              r_zy = c(1, 1))
  r_TstarU <- 0
  k <- 1
  expected_answer <- c(sqrt(10), sqrt(10))
  ans <- get_M(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_M properly computes M for multiple sims each with corresponding
           values of r_TstarU and k", {
  obs <- list(r_Ty = c(sqrt(3/4), sqrt(3/4)),
              r_Tz = c(0, 0),
              r_zy = c(1, 1))
  r_TstarU <- c(0, 0)
  k <- c(1, 1)
  expected_answer <- c(sqrt(10), sqrt(10))
  ans <- get_M(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_M properly computes M for one sim across a range of
          values of r_TstarU and k", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 0,
              r_zy = 1)
  r_TstarU <- c(0, 0)
  k <- c(1, 1)
  expected_answer <- c(sqrt(10), sqrt(10))
  ans <- get_M(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_M throws an error for mismatched dimensions", {
  obs <- list(r_Ty = c(sqrt(3/4), sqrt(3/4)),
              r_Tz = c(0, 0),
              r_zy = c(1, 1))
  r_TstarU <- c(0, 0, 0)
  k <- c(1, 1, 1)
  expect_error(get_M(r_TstarU, k, obs))
})

test_that("get_s_u properly computes s_u", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 0,
              r_zy = 1,
              s2_y = 1)
  r_TstarU <- 0
  k <- 1
  expected_answer <- 1 / 2
  ans <- get_s_u(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_s_u properly computes s_u for multiple simulations", {
  obs <- list(r_Ty = c(sqrt(3/4), sqrt(3/4), sqrt(3/4)),
              r_Tz = c(0, 0, 0),
              r_zy = c(1, 1, 1),
              s2_y = c(1, 1, 1))
  r_TstarU <- 0
  k <- 1
  expected_answer <- c(1/2, 1/2, 1/2)
  ans <- get_s_u(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_s_u properly computes s_u for multiple simulations and
          corresponding values of r_TstarU and k", {
  obs <- list(r_Ty = c(sqrt(3/4), sqrt(3/4), sqrt(3/4)),
              r_Tz = c(0, 0, 0),
              r_zy = c(1, 1, 1),
              s2_y = c(1, 1, 1))
  r_TstarU <- c(0, 0, 0)
  k <- c(1, 1, 1)
  expected_answer <- c(1/2, 1/2, 1/2)
  ans <- get_s_u(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_s_u properly computes s_u for one simulation and across a grid of
          values of r_TstarU and k", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 0,
              r_zy = 1,
              s2_y = 1)
  r_TstarU <- c(0, 0, 0)
  k <- c(1, 1, 1)
  expected_answer <- c(1/2, 1/2, 1/2)
  ans <- get_s_u(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_s_u throws an error if multiple simulations do not have
          corresponding values of r_TstarU and k", {
  obs <- list(r_Ty = rep(sqrt(3/4), 3),
              r_Tz = rep(0, 3),
              r_zy = rep(1, 3),
              s2_y = rep(1, 3))
  r_TstarU <- rep(0, 2)
  k <- rep(1, 2)
  expect_error(get_s_u(r_TstarU, k, obs))
})

test_that("get_beta properly computes beta", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 1,
              r_zy = 1,
              s2_y = 1,
              s2_T = 1)
  r_TstarU <- 0
  k <- 1
  expected_answer <- sqrt(3) / 2
  ans <- get_beta(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta properly computes beta for multiple simulations", {
  obs <- list(r_Ty = rep(sqrt(3/4), 3),
              r_Tz = rep(1, 3),
              r_zy = rep(1, 3),
              s2_y = rep(1, 3),
              s2_T = rep(1, 3))
  r_TstarU <- 0
  k <- 1
  expected_answer <- rep(sqrt(3) / 2, 3)
  ans <- get_beta(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta properly computes beta for multiple simulations each
          with a corresponding r_TstarU and k.", {
  obs <- list(r_Ty = rep(sqrt(3/4), 3),
              r_Tz = rep(1, 3),
              r_zy = rep(1, 3),
              s2_y = rep(1, 3),
              s2_T = rep(1, 3))
  r_TstarU <- rep(0, 3)
  k <- rep(1, 3)
  expected_answer <- rep(sqrt(3) / 2, 3)
  ans <- get_beta(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta properly computes beta for one simulation across a grid of
          r_TstarU and k values.", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 1,
              r_zy = 1,
              s2_y = 1,
              s2_T = 1)
  r_TstarU <- rep(0, 3)
  k <- rep(1, 3)
  expected_answer <- rep(sqrt(3) / 2, 3)
  ans <- get_beta(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta throws an error when multiple simulations do not have
          a corresponding r_TstarU and k.", {
  obs <- list(r_Ty = rep(sqrt(3/4), 3),
              r_Tz = rep(1, 3),
              r_zy = rep(1, 3),
              s2_y = rep(1, 3),
              s2_T = rep(1, 3))
  r_TstarU <- rep(0, 2)
  k <- rep(1, 2)
  expect_error(get_beta(r_TstarU, k, obs))
})

test_that("get_beta_lower properly computes bound", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 1,
              r_zy = 1,
              s2_y = 1,
              s2_T = 1)
  r_TstarU_max <- 0
  k_min <- sqrt(3/4)
  k_max <- 1
  expected_answer <- sqrt(3/4)
  ans <- get_beta_lower(r_TstarU_max, k_min, k_max, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta_lower properly computes bound for multiple simulations but
          same restrictions", {
  obs <- list(r_Ty = rep(sqrt(3/4), 3),
              r_Tz = rep(1, 3),
              r_zy = rep(1, 3),
              s2_y = rep(1, 3),
              s2_T = rep(1, 3))
  r_TstarU_max <- 0
  k_min <- sqrt(3/4)
  k_max <- 1
  expected_answer <- rep(sqrt(3/4), 3)
  ans <- get_beta_lower(r_TstarU_max, k_min, k_max, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta_lower properly computes bound for multiple simulations each
          with corresponding restrictions", {
  obs <- list(r_Ty = rep(sqrt(3/4), 3),
              r_Tz = rep(1, 3),
              r_zy = rep(1, 3),
              s2_y = rep(1, 3),
              s2_T = rep(1, 3))
  r_TstarU_max <- rep(0, 3)
  k_min <- rep(sqrt(3/4), 3)
  k_max <- rep(1, 3)
  expected_answer <- rep(sqrt(3/4), 3)
  ans <- get_beta_lower(r_TstarU_max, k_min, k_max, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta_lower properly computes bound for one simulation across a
          grid of restrictions", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 1,
              r_zy = 1,
              s2_y = 1,
              s2_T = 1)
  r_TstarU_max <- rep(0, 3)
  k_min <- rep(sqrt(3/4), 3)
  k_max <- rep(1, 3)
  expected_answer <- rep(sqrt(3/4), 3)
  ans <- get_beta_lower(r_TstarU_max, k_min, k_max, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta_lower throws an error dimensions mismatch", {
  obs <- list(r_Ty = c(sqrt(3/4), sqrt(3/4)),
              r_Tz = c(1, 1),
              r_zy = c(1, 1),
              s2_y = c(1, 1),
              s2_T = c(1, 1))
  r_TstarU_max <- rep(0, 3)
  k_min <- rep(sqrt(3/4), 3)
  k_max <- rep(1, 3)
  expect_error(get_beta_lower(r_TstarU_max, k_min, k_max, obs))
})

test_that("get_beta_upper properly computes bound", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 1,
              r_zy = 1,
              s2_y = 1,
              s2_T = 1)
  r_TstarU_max <- 0
  k_min <- sqrt(3/4)
  k_max <- 1
  expected_answer <- 1
  ans <- get_beta_upper(r_TstarU_max, k_min, k_max, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta_upper properly computes bound for multiple simulations but
          same restrictions", {
  obs <- list(r_Ty = rep(sqrt(3/4), 3),
              r_Tz = rep(1, 3),
              r_zy = rep(1, 3),
              s2_y = rep(1, 3),
              s2_T = rep(1, 3))
  r_TstarU_max <- 0
  k_min <- sqrt(3/4)
  k_max <- 1
  expected_answer <- rep(1, 3)
  ans <- get_beta_upper(r_TstarU_max, k_min, k_max, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta_upper properly computes bound for multiple simulations each
          with corresponding restrictions", {
  obs <- list(r_Ty = rep(sqrt(3/4), 3),
              r_Tz = rep(1, 3),
              r_zy = rep(1, 3),
              s2_y = rep(1, 3),
              s2_T = rep(1, 3))
  r_TstarU_max <- rep(0, 3)
  k_min <- rep(sqrt(3/4), 3)
  k_max <- rep(1, 3)
  expected_answer <- rep(1, 3)
  ans <- get_beta_upper(r_TstarU_max, k_min, k_max, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta_upper properly computes bound for one simulation across a
          grid of restrictions", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 1,
              r_zy = 1,
              s2_y = 1,
              s2_T = 1)
  r_TstarU_max <- rep(0, 3)
  k_min <- rep(sqrt(3/4), 3)
  k_max <- rep(1, 3)
  expected_answer <- rep(1, 3)
  ans <- get_beta_upper(r_TstarU_max, k_min, k_max, obs)
  expect_equal(expected_answer, ans)
})

test_that("get_beta_upper throws an error dimensions mismatch", {
  obs <- list(r_Ty = c(sqrt(3/4), sqrt(3/4)),
              r_Tz = c(1, 1),
              r_zy = c(1, 1),
              s2_y = c(1, 1),
              s2_T = c(1, 1))
  r_TstarU_max <- rep(0, 3)
  k_min <- rep(sqrt(3/4), 3)
  k_max <- rep(1, 3)
  expect_error(get_beta_upper(r_TstarU_max, k_min, k_max, obs))
})
