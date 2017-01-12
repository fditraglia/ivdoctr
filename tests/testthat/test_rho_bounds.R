context("Rho bounds")

test_that("candidate1 computes properly 1", {
  obs <- list(r_Ty = 0,
              r_Tz = 0,
              r_zy = 0)
  r_TstarU_upper <- 1
  r_TstarU_lower <- 0
  k_upper <- 1
  k_lower <- 0.1
  expected_answer <- list(min_corner = 0, max_corner = 0)
  ans <- candidate1(r_TstarU_lower, r_TstarU_upper, k_lower, k_upper, obs)
  expect_equal(expected_answer, ans)
})

test_that("candidate1 is appropriately vectorized wrt obs", {
  obs <- list(r_Ty = c(0, 0),
              r_Tz = c(0, 0),
              r_zy = c(0, 0))
  r_TstarU_upper <- 1
  r_TstarU_lower <- 0
  k_upper <- 1
  k_lower <- 0.1
  expected_answer <- list(min_corner = c(0, 0), max_corner = c(0, 0))
  ans <- candidate1(r_TstarU_lower, r_TstarU_upper, k_lower, k_upper, obs)
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
  expected_answer <- list(min_edge = 1, max_edge = 1)
  ans <- candidate2(r_TstarU_lower, r_TstarU_upper, k_lower, k_upper, obs)
  expect_equivalent(expected_answer, ans)
})

test_that("candidate2 is vectorized wrt obs", {
  obs <- list(r_Ty = c(0, 0, 0),
              r_Tz = c(0, 0, 0),
              r_zy = c(1, 1, 1))
  r_TstarU_upper <- 1
  r_TstarU_lower <- 0
  k_upper <- 1
  k_lower <- 0.1
  expected_answer <- list(min_edge = c(1, 1, 1), max_edge = c(1, 1, 1))
  ans <- candidate2(r_TstarU_lower, r_TstarU_upper, k_lower, k_upper, obs)
  expect_equivalent(expected_answer, ans)
})

test_that("candidate3 properly computes zeros of polynomial to get kappas", {
  k_lower <- 0.8
  k_upper <- 0.9
  r_TstarU_lower <- -0.8
  r_TstarU_upper <- -0.4
  obs <- list(r_Ty = -0.9,
              r_Tz = 0.24,
              r_zy = -0.3)

  ans <- candidate3(r_TstarU_lower, r_TstarU_upper, k_lower, k_upper, obs)
  test1 <- with(obs, -r_TstarU_upper ^ 2 * r_Tz ^ 2 *
                      (ans$k_roots$real_roots_upper - r_Ty ^ 2) ^ 3 +
                     ((2 * r_Ty * r_Tz - r_zy * r_Ty ^ 2) *
                      ans$k_roots$real_roots_upper -
                      r_Tz * r_Ty ^ 3) ^ 2 * (1 - r_TstarU_upper ^ 2))
  test2 <- with(obs, -r_TstarU_lower ^ 2 * r_Tz ^ 2 *
                     (ans$k_roots$real_roots_lower - r_Ty ^ 2) ^ 3 +
                     ((2 * r_Ty * r_Tz - r_zy * r_Ty ^ 2) *
                     ans$k_roots$real_roots_lower -
                     r_Tz * r_Ty ^ 3) ^ 2 * (1 - r_TstarU_lower ^ 2))
  expect_equal(sum(test1, na.rm = TRUE), 0, tol = 10 ^ -4)
  expect_equal(sum(test2, na.rm = TRUE), 0, tol = 10 ^ -4)
})

test_that("candidate3 properly computes zeros of polynomial to get kappas with
          multiple simulations of Sigma", {
  k_lower <- rep(0.8, 4)
  k_upper <- rep(0.9, 4)
  r_TstarU_lower <- rep(-0.8, 4)
  r_TstarU_upper <- rep(-0.4, 4)
  obs <- list(r_Ty = c(-0.9, -0.9, -0.9, -0.9),
              r_Tz = c(0.24, 0.24, 0.24, 0.24),
              r_zy = c(-0.3, -0.3, -0.3, -0.3))

  ans <- candidate3(r_TstarU_lower, r_TstarU_upper, k_lower, k_upper, obs)
  test1 <- with(obs, -r_TstarU_upper ^ 2 * r_Tz ^ 2 *
                      (ans$k_roots$real_roots_upper - r_Ty ^ 2) ^ 3 +
                     ((2 * r_Ty * r_Tz - r_zy * r_Ty ^ 2) *
                      ans$k_roots$real_roots_upper -
                      r_Tz * r_Ty ^ 3) ^ 2 * (1 - r_TstarU_upper ^ 2))
  test2 <- with(obs, -r_TstarU_lower ^ 2 * r_Tz ^ 2 *
                      (ans$k_roots$real_roots_lower - r_Ty ^ 2) ^ 3 +
                     ((2 * r_Ty * r_Tz - r_zy * r_Ty ^ 2) *
                      ans$k_roots$real_roots_lower -
                      r_Tz * r_Ty ^ 3) ^ 2 * (1 - r_TstarU_lower ^ 2))
  expect_equal(sum(test1, na.rm = TRUE), 0, tol = 10 ^ -4)
  expect_equal(sum(test2, na.rm = TRUE), 0, tol = 10 ^ -4)
})

# test_that("alternative full solution", {
#   set.seed(1234)
#   nsim <- 500
#   lengthGrid <- 500
#
#   # Generating a bunch of correlation draws
#   sims_rho_Tz <- 2 * runif(nsim) - 1
#   sims_rho_Ty <- 2 * runif(nsim) - 1
#   sims_rho_zy <- 2 * runif(nsim) - 1
#   obs <- list(r_Ty = sims_rho_Ty,
#               r_zy = sims_rho_zy,
#               r_Tz = sims_rho_Tz,
#               T_Rsq = rep(0, nsim))
#
#   # Getting unrestricted bounds for kappa
#   k_bounds <- get_k_bounds_unrest(obs, tilde = FALSE)
#
#   # Generating matrix that compares numeric and analytic solutions
#   results <- matrix(NA, nsim, 4)
#   colnames(results) <- c("minNumeric", "minAnalytic", "maxNumeric", "maxAnalytic")
#
#   # Generating matrix that records the relevant arguments for the solutions
#   args_results <- matrix(NA, nsim, 10)
#   colnames(args_results) <- c("K_L", "K_U", "K_MIN", "K_MAX",
#                               "r_TstarU_L", "r_TstarU_U", "r_TstarU_MIN", "r_TstarU_MAX",
#                               "Set solution MIN", "Set solution MAX")
#
#   test_positivedefiniteinputs <- NULL
#
#   ### Loop over every simulation
#
#   for (i in 1:nsim) {
#
#     Klower <- k_bounds$Lower[i]
#     Kupper <- k_bounds$Upper[i]
#     sim_obs <- list(r_Ty = sims_rho_Ty[i],
#                     r_Tz = sims_rho_Tz[i],
#                     r_zy = sims_rho_zy[i])
#     Sigma <- matrix(c(1, sim_obs$r_Ty, sim_obs$r_Tz,
#                       sim_obs$r_Ty, 1, sim_obs$r_zy,
#                       sim_obs$r_Tz, sim_obs$r_zy, 1), nrow = 3, ncol = 3)
#     test_positivedefiniteinputs[i] <- matrixcalc::is.positive.definite(Sigma)
#
#     # Generating a bunch of draws for user-imposed bounds on kappa.
#     K_1 <- (1 - Klower) * runif(1) + Klower
#     K_2 <- (1 - Klower) * runif(1) + Klower
#     r_TstarU_1 <- 2 * runif(1) - 1
#     r_TstarU_2 <- 2 * runif(1) - 1
#
#     # Creating the user-imposed upper and lower bounds on k and r_TstarU
#     K_L <- min(K_1, K_2)
#     K_U <- max(K_1, K_2)
#     K_bounds <- c(K_L, K_U)
#     r_TstarU_L <- min(r_TstarU_1, r_TstarU_2)
#     r_TstarU_U <- max(r_TstarU_1, r_TstarU_2)
#     r_TstarU_bounds <- c(r_TstarU_L, r_TstarU_U)
#
#     # If covariance matrix is PSD and user lower bound is binding, then we
#     # compute the solution
#     if (test_positivedefiniteinputs[i] & K_L >= Klower) {
#
#       # Create grid of kappa and r_TstarU and then evaluate to get r_uz
#       K <- seq(K_L, K_U, length.out = lengthGrid)
#       r_TstarU <- seq(r_TstarU_L, r_TstarU_U, length.out = lengthGrid)
#       Rzumatrix <- matrix(, nrow = lengthGrid, ncol = lengthGrid)
#       for (i in 1:lengthGrid) {
#         for (j in 1:lengthGrid) {
#           Rzumatrix[i, j] <- get_r_uz(r_TstarU[i], K[j], sim_obs)
#         }
#       }
#       Rzumatrix <- outer(r_TstarU, K, function(r_TstarU, K) get_r_uz(r_TstarU, K, sim_obs))
#
#       # Store max and min r_zu
#       results[i, "minNumeric"] <- min(Rzumatrix)
#       results[i, "maxNumeric"] <- max(Rzumatrix)
#
#       # Store arguments that generated those max and min values as well as the bounds
#       args_results[i, "K_L"] <- K_L
#       args_results[i, "K_U"] <- K_U
#       args_results[i, "K_MIN"] <- K[which(Rzumatrix == min(Rzumatrix),
#                                           arr.ind = TRUE)[2]]
#       args_results[i, "K_MAX"] <- K[which(Rzumatrix == max(Rzumatrix),
#                                           arr.ind = TRUE)[2]]
#       args_results[i, "r_TstarU_L"] <- r_TstarU_L
#       args_results[i, "r_TstarU_U"] <- r_TstarU_U
#       args_results[i, "r_TstarU_MIN"] <- r_TstarU[which(Rzumatrix == min(Rzumatrix),
#                                                         arr.ind = TRUE)[1]]
#       args_results[i, "r_TstarU_MAX"] <- r_TstarU[which(Rzumatrix == max(Rzumatrix),
#                                                         arr.ind = TRUE)[1]]
#
#       ## MIN RESULTS
#       # Corners
#       if (args_results[i, "K_MIN"] %in% c(K_L, K_U) &
#          args_results[i, "r_TstarU_MIN"] %in% c(r_TstarU_L, r_TstarU_U)) {
#         args_results[i, "Set solution MIN"] <- "Corner"
#
#       # Interior kappa
#       } else if (args_results[i, "r_TstarU_MIN"] %in% c(r_TstarU_L, r_TstarU_U)) {
#         args_results[i, "Set solution MIN"] <- "Left or right edge"
#
#       # Interior r_TstarU
#       } else if (args_results[i, "K_MIN"] %in% c(K_L, K_U)) {
#         args_results[i, "Set solution MIN"] <- "Top or bottom edge"
#
#       # Interior
#       } else {
#         args_results[i, "Set solution MIN"] <- "Interior"
#       }
#
#       ## MAX RESULTS
#       # Corners
#       if (args_results[i, "K_MAX"] %in% c(K_L, K_U) &
#          args_results[i, "r_TstarU_MAX"] %in% c(r_TstarU_U, r_TstarU_L)) {
#         args_results[i, "Set solution MAX"] <- "Corner"
#
#       # Interior kappa
#       } else if (args_results[i, "r_TstarU_MAX"] %in% c(r_TstarU_L, r_TstarU_U)) {
#         args_results[i, "Set solution MAX"] <- "Left or right edge"
#
#       # Interior r_TstarU
#       } else if (args_results[i, "K_MAX"] %in% c(K_L, K_U)) {
#         args_results[i, "Set solution MAX"] <- "Top or bottom edge"
#
#       # Interior
#       } else {
#         args_results[i, "Set solution MAX"] <- "Interior"
#       }
#
#       # Calculate roots using polynomial
#       d <- with(sim_obs, -r_Tz ^ 2 * r_TstarU_bounds ^ 2)
#       c <- with(sim_obs, 3 * r_Ty ^ 2 * (r_Tz ^ 2) * (r_TstarU_bounds ^ 2) +
#                          ((2 * r_Ty * r_Tz - r_zy * (r_Ty ^ 2)) ^ 2) *
#                          (1 - r_TstarU_bounds ^ 2))
#       b <- with(sim_obs, -3 * (r_Ty ^ 4) * (r_Tz ^ 2) * (r_TstarU_bounds ^ 2) -
#                           2 * (2 * r_Ty * r_Tz - r_zy * (r_Ty ^ 2)) *
#                           (r_Tz * (r_Ty ^ 3)) * (1 - r_TstarU_bounds ^ 2))
#       a <- rep(with(sim_obs, (r_Ty ^ 6) * (r_Tz ^ 2)), 2)
#
#       all_roots <- apply(cbind(a, b, c, d), 1, polyroot)
#       roots <- Re(all_roots)
#       roots_L <- roots[, 1]
#       roots_L <- roots_L[(roots_L >= K_L) & (roots_L <= K_U)]
#       roots_U <- roots[, 2]
#       roots_U <- roots_U[(roots_U >= K_L) & (roots_U <= K_U)]
#
#       S1_L <- matrix(NA, length(roots_L), 3)
#       S1_U <- matrix(NA, length(roots_U), 3)
#
#       S1_L[, 1] <- roots_L
#       S1_L[, 2] <- r_TstarU_L
#
#       S1_U[, 1] <- roots_U
#       S1_U[, 2] <- r_TstarU_U
#
#       ### CALCULATE THE SET S2; Corner solutions
#       S2_corners <- matrix(NA, 4, 3)
#       S2_corners[1, ] <- c(K_L, r_TstarU_L, NA)
#       S2_corners[2, ] <- c(K_L, r_TstarU_U, NA)
#       S2_corners[3, ] <- c(K_U, r_TstarU_L, NA)
#       S2_corners[4, ] <- c(K_U, r_TstarU_U, NA)
#
#       ### CALCULATE THE SET S3: Edges with extreme points of K
#       a_L <- with(sim_obs, -r_Tz * sqrt(K_L - (r_Ty ^ 2)) / (r_Ty * r_Tz - K_L * r_zy))
#       a_U <- with(sim_obs, -r_Tz * sqrt(K_U - (r_Ty ^ 2)) / (r_Ty * r_Tz - K_U * r_zy))
#
#       S3_int <- matrix(NA, 2, 3)
#
#       if ((a_L / sqrt(1 + a_L ^ 2) >= r_TstarU_L) &
#           (a_L / sqrt(1 + a_L ^ 2) <= r_TstarU_U)) {
#         S3_int[1, ] <- c(K_L, a_L / sqrt(1 + a_L ^ 2), NA)
#       }
#
#       if ((a_U / sqrt(1 + a_U ^ 2) >= r_TstarU_L) &
#           (a_U / sqrt(1 + a_U ^ 2) <= r_TstarU_U)) {
#         S3_int[2, ] <- c(K_U, a_U / sqrt(1 + a_U ^ 2),NA)
#       }
#
#       ### COMPUTE CANDIDATE SET
#       S <- rbind(S1_L, S1_U, S2_corners, S3_int)
#
#       for (m in (1:nrow(S))) {
#         if (!is.na(S[m, 1]) & !is.na(S[m, 2])) {
#           S[m, 3] <- get_r_uz(S[m, 2], S[m, 1], sim_obs)
#         }
#       }
#
#       ## Output results
#       results[i, "minAnalytic"] <- min(S[, 3], na.rm = TRUE)
#       results[i, "maxAnalytic"] <- max(S[, 3], na.rm = TRUE)
#
#       if (results[i, "minAnalytic"] > results[i, "minNumeric"]) {
#         print("min not working")
#       }
#
#       if (results[i, "maxAnalytic"] < results[i, "maxNumeric"]) {
#         print("max not working")
#       }
#     }
#   }
#
#   ### Figure out whether the analytical minimum is less than the lowest value in the grid
#   ### and whether the max is larger than the value in the grid.
#
#   diff_min <- results[, "minNumeric"] - results[, "minAnalytic"]
#   mean(diff_min >= 0, na.rm = TRUE)
#   expect_equal(0, length(args_results[diff_min < 0 & !is.na(diff_min), ]))
#
#   diff_max <- results[, "maxNumeric"] - results[, "maxAnalytic"]
#   mean(diff_max <= 0, na.rm = TRUE)
#   expect_equal(0, length(args_results[diff_max > 0 & !is.na(diff_max), ]))
# })
