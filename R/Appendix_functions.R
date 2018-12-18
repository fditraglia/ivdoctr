#' G function from Proposition A.2
#'
#' @param kappa Kappa value
#' @param r_TstarU r_TstarU value
#' @param obs_draws a row of the data.frame of observable draws
#'
#' @return G value
#' @export
#'
g_functionA2 <- function(kappa, r_TstarU, obs_draws) {
  t1 <- with(obs_draws, sqrt(s2_y) / kappa)
  t2 <- with(obs_draws, r_Tz * sqrt(kappa - r_Ty ^ 2) * r_TstarU / (sqrt(1 - r_TstarU ^ 2)))
  t3 <- with(obs_draws, r_Ty * r_Tz - kappa * r_zy)
  ans <- t1 * (t2 - t3)
  return(ans)
}

#' B function from Proposition A3
#'
#' @param kappa Kappa value
#' @param obs_draws Row of the data.frame of observable draws
#' @param p Treatment probability from data
#'
#' @return A min and a max of the B function
#' @export
#'
b_functionA3 <- function(kappa, obs_draws, p) {
  g_bounds <- g_functionA2(kappa, c(-1, 1), obs_draws)
  psi_bounds <- c(get_psi_lower(obs_draws$s2_T, p, kappa),
                  get_psi_upper(obs_draws$s2_T, p, kappa))
  B_args <- expand.grid(g = g_bounds, psi = psi_bounds)
  B_vals <- with(B_args, (1 + psi) * (s_zy / s_Tz - g))
  ans <- c(min(B_vals), max(B_vals))
  return(ans)
}

#' Returns beta bounds in binary case using grid search
#'
#' @param obs_draws Row of the data.frame of observable draws
#' @param p Treatment probability from data
#'
#' @return Min and max values for beta
#' @export
#'
get_beta_bounds_binary <- function(obs_draws, p) {
  L <- get_L(obs_draws)
  for (i in 1:length(L)) {
    kappas <- matrix(seq(L[i], 1, by = 0.0001), ncol = 1)
    vals <- apply(kappas, 1, b_functionA3, obs_draws = obs_draws, p = p)

  }
  ans <- c(min(vals), max(vals))
  return(ans)
}
