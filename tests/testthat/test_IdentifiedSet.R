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

test_that("get_r_uz is vectorized wrt r_TstarU and k", {
  obs <- list(r_Ty = sqrt(3/4),
              r_Tz = 0,
              r_zy = 1)
  r_TstarU <- c(0, 0)
  k <- c(1, 1)
  expected_answer <- c(2, 2)
  ans <- get_r_uz(r_TstarU, k, obs)
  expect_equal(expected_answer, ans)
})

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
