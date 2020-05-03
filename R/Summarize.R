#' @importFrom AER ivreg
#' @importFrom coda HPDinterval as.mcmc
#' @importFrom sandwich vcovHC
NULL

#' Computes OLS and IV estimates
#' @param y_name Character vector of the name of the dependent variable
#' @param T_name Character vector of the names of the preferred regressors
#' @param z_name Character vector of the names of the instrumental variables
#' @param data Data to be analyzed
#' @param controls Character vector containing the names of the exogenous regressors
#' @param robust Boolean of whether to compute heteroskedasticity-robust standard errors
#' @return List of beta estimates and associated standard errors for OLS and
#'   IV estimation
get_estimates <- function(y_name, T_name, z_name, data, controls = NULL,
                          robust = FALSE) {
  first_stage <- reformulate(c(z_name, controls), response = NULL)
  second_stage <- reformulate(c(T_name, controls), response = y_name)
  OLS <- lm(second_stage, data)
  IV <-  AER::ivreg(second_stage, first_stage, data)

  b_OLS <- coef(OLS)[T_name]
  b_IV <- coef(IV)[T_name]

  if (robust) {
    se_OLS <- sqrt(diag(sandwich::vcovHC(OLS, type = 'HC0')))[T_name]
    se_IV <- sqrt(diag(sandwich::vcovHC(IV, type = 'HC0')))[T_name]
  } else {
    se_OLS <- sqrt(diag(vcov(OLS)))[T_name]
    se_IV <- sqrt(diag(vcov(IV)))[T_name]
  }

  list(n = nrow(data),
       b_OLS = b_OLS,
       se_OLS = se_OLS,
       b_IV = b_IV,
       se_IV = se_IV)
}

get_HPDI <- function(draws, level = 0.9) {
  interval <- coda::HPDinterval(coda::as.mcmc(draws), level)
  lower <- interval[[1]]
  upper <- interval[[2]]
  return(data.frame(lower = lower, median = median(draws, na.rm = TRUE), upper = upper))
}

summarize_bounds <- function(draws) {
  unrestricted <- with(draws$unrestricted,
                       rbind(k_lower = get_HPDI(k_lower),
                             r_uz_lower = get_HPDI(r_uz_lower),
                             r_uz_upper = get_HPDI(r_uz_upper)))
  p_valid <- get_p_valid(draws)
  list(unrestricted = unrestricted,
       r_TstarU_restriction = draws$r_TstarU_restriction,
       k_restriction = draws$k_restriction,
       p_empty = mean(draws$empty),
       p_valid = p_valid)
}

summarize_posterior <- function(draws) {
  HPDI <- with(draws$posterior, rbind(r_uz = get_HPDI(r_uz),
                                      beta = get_HPDI(beta)))
  list(r_TstarU_restriction = draws$r_TstarU_restriction,
       k_restriction = draws$k_restriction,
       HPDI = HPDI)
}

summarize_posterior_binary <- function(draws, p) {
  # Calculating new beta with psi adjustment
  kappa <- draws$posterior$k
  s2_T <- draws$observables$s2_T
  psi_lower <- get_psi_lower(s2_T, p, kappa)
  psi_upper <- get_psi_upper(s2_T, p, kappa)
  psi_bounds <- cbind(psi_lower, psi_upper)
  psi_draw <- apply(psi_bounds, 1, function(x) runif(1, x[1], x[2]))
  new_beta <- draws$posterior$beta * (1 + psi_draw)
  HPDI <- with(draws$posterior, rbind(r_uz = get_HPDI(r_uz),
                                      beta = get_HPDI(new_beta)))
  list(r_TstarU_restriction = draws$r_TstarU_restriction,
       k_restriction = draws$k_restriction,
       HPDI = HPDI)
}
