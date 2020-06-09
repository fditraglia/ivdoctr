#' Computes beliefs that support valid instrument
#'
#' @param post_draws data.frame of posterior draws
#' @param obs_draws data.frame of draws of reduced form parameters
#'
#' @return data.frame of new draws
#'
get_new_draws <- function(obs_draws, post_draws) {
  kappa <- post_draws$k
  A <- get_A(obs_draws, kappa)
  r_TstarU <- sqrt(A ^ 2 / (1 + A ^ 2)) * sign(A)
  new_draws <- data.frame(r_TstarU = r_TstarU,
                          k = kappa,
                          r_uz = rep(0, length(kappa)),
                          s_u = post_draws$s_u,
                          post_draws = post_draws$beta)
  return(new_draws)
}

get_A <- function(obs_draws, kappa) {
  num <- with(obs_draws, r_Ty * r_Tz - r_zy * kappa)
  denom <- with(obs_draws, r_Tz * sqrt(kappa - r_Ty ^ 2))
  ans <- num / denom
  return(ans)
}
