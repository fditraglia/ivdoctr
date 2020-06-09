#' Computes the upper bound of psi for binary data
#'
#' @param s2_T Vector of s2_T draws from observables
#' @param p Treatment probability from binary data
#' @param kappa Vector of kappa, NOTE: kappa_tilde in the paper
#'
#' @return Vector of upper bounds for psi
#'
get_psi_upper <- function(s2_T, p, kappa) {
  ans <- (-s2_T * (1 - kappa)) / max(p, 1 - p)
  return(ans)
}

#' Computes the lower bound of psi for binary data
#'
#' @param s2_T Vector of s2_T draws from observables
#' @param p Treatment probability from binary data
#' @param kappa Vector of kappa, NOTE: kappa_tilde in the paper
#'
#' @return Vector of lower bounds for psi
#'
get_psi_lower <- function(s2_T, p, kappa) {
  m <- max((1 - p) * (2 * p - 2), p * (1 - 2 * p))
  ind <- (s2_T * (1 - kappa) <= m)
  psi1 <- (-s2_T * (1 - kappa)) / min(p, 1 - p)
  psi2 <- 2 * sqrt(p * (1 - p) - s2_T * (1 - kappa)) - 1
  psi_lower <- psi1 * ind + psi2 * (1 - ind)
  return(psi_lower)
}

#' Computes L, lower bound for kappa_tilde in paper
#'
#' @param draws data.frame of observables of simulated data
#'
#' @return Vector of L values
#'
get_L <- function(draws) {
  num <- with(draws, r_Ty ^ 2 + r_Tz ^ 2 - 2 * r_Ty * r_Tz * r_zy)
  denom <- with(draws, 1 - r_zy ^ 2)
  L <- num / denom
  return(L)
}

#' Computes a0 and a1 bounds
#'
#' @param draws data.frame of observables of simulated data
#' @param p Treatment probability from binary data
#'
#' @return List of alpha bounds
#'
get_alpha_bounds <- function(draws, p) {
  L <- get_L(draws)
  a0 <- draws$s2_T * (1 - L) / (1 - p)
  a1 <- draws$s2_T * (1 - L) / p
  ans <- list(a0 = a0, a1 = a1)
  return(ans)
}
