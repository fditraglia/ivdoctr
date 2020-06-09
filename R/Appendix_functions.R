#' G function from Proposition A.2
#'
#' @param kappa Kappa value
#' @param r_TstarU r_TstarU value
#' @param obs_draws a row of the data.frame of observable draws
#'
#' @return G value
#'
g_functionA2 <- function(kappa, r_TstarU, obs_draws) {
  t1 <- with(obs_draws, sqrt(s2_z) / s_Tz * sqrt(s2_y) / kappa)
  t2 <- with(obs_draws, r_Tz * sqrt(kappa - r_Ty ^ 2) * r_TstarU / (sqrt(1 - r_TstarU ^ 2)))
  t3 <- with(obs_draws, r_Ty * r_Tz - kappa * r_zy)
  ans <- t1 * (t2 - t3)
  return(ans)
}

#' B function from Proposition A3
#'
#' @param obs_draws Row of the data.frame of observable draws
#' @param g Value from g function
#' @param psi Psi value
#'
#' @return A min and a max of the B function
#'
b_functionA3 <- function(obs_draws, g, psi) {
  vals <- with(obs_draws, (1 + psi) * (s_zy / s_Tz - g))
  ans <- cbind(beta_lower = min(vals), beta_upper = max(vals))
  return(ans)
}

#' Returns beta bounds in binary case using grid search
#'
#' @param obs_draws Row of the data.frame of observable draws
#' @param p Treatment probability from data
#' @param r_TstarU_restriction 2-element vector of restrictions on r_TstarU
#'
#' @return Min and max values for beta
#'
get_beta_bounds_binary <- function(obs_draws, p, r_TstarU_restriction) {
  L <- get_L(obs_draws)
  beta_bounds <- NULL
  for (i in 1:length(L)) {
    kappas <- matrix(seq(L[i], 1, length.out = 100), ncol = 1)
    for (j in 1:length(kappas)) {
      psi_lower <- get_psi_lower(obs_draws$s2_T[i], p, kappas[j])
      psi_upper <- get_psi_upper(obs_draws$s2_T[i], p, kappas[j])
      psi <- c(psi_lower, psi_upper)
      g_1 <- g_functionA2(kappas[j], r_TstarU_restriction[1], obs_draws[i, ])
      g_2 <- g_functionA2(kappas[j], r_TstarU_restriction[2], obs_draws[i, ])
      g_val <- c(min(g_1, g_2), max(g_1, g_2))
      inputs <- expand.grid(psi = psi, g = g_val)
      beta_bounds <- rbind(beta_bounds,
                           b_functionA3(obs_draws = obs_draws[i, ], g = inputs$g, psi = inputs$psi))
    }
  }
  return(beta_bounds)
}

#' Generates beta bounds off of beta draws
#'
#' @param draws Posterior draws
#' @param n_observables Number of observable draws
#'
#' @return Upper and lower bounds of beta based on posterior draws
#'
get_beta_bounds_binary_post <- function(draws, n_observables) {
  reps <- nrow(draws$posterior) / nrow(draws$observables)
  beta_bounds <- NULL
    for (i in 1:n_observables) {
      start <- (i - 1) * reps + 1
      end <- i * reps
      beta_lower <- min(draws$posterior$beta[start:end])
      beta_upper <- max(draws$posterior$beta[start:end])
      beta_bounds <- rbind(beta_bounds, c(beta_lower, beta_upper))
    }
  return(beta_bounds)
}
