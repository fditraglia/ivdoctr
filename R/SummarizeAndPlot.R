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
  upper<- interval[[2]]
  return(data.frame(lower = lower, median = median(draws), upper = upper))
}

summarize_bounds <- function(draws) {
  unrestricted <- with(draws$unrestricted,
                       rbind(k_lower = get_HPDI(k_lower),
                             r_uz_lower = get_HPDI(r_uz_lower),
                             r_uz_upper = get_HPDI(r_uz_upper)))
  # Add in restricted r_uz_lower and upper bounds later:
  restricted <- with(draws$restricted,
                     rbind(beta_lower = get_HPDI(beta_lower),
                           beta_upper = get_HPDI(beta_upper)))
  list(unrestricted = unrestricted,
       r_TstarU_restriction = draws$r_TstarU_restriction,
       k_restriction = draws$k_restriction,
       p_empty = mean(draws$empty),
       restricted = restricted)
}

summarize_posterior <- function(draws) {
  HPDI <- with(draws$posterior, rbind(r_uz = get_HPDI(r_uz),
                                      beta = get_HPDI(beta)))
  list(r_TstarU_restriction = draws$r_TstarU_restriction,
       k_restriction = draws$k_restriction,
       p_empty = mean(draws$empty),
       HPDI = HPDI)
}



