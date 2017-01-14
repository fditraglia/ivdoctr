# #----------------------- Colonial Origins -------------------------------------
# load("./data/colonial.rda")
#
# # Part I of table
# summary_1 <- get_estimates("logpgp95", "avexpr", "logem4", colonial)
# obs <- get_observables("logpgp95", "avexpr", "logem4", colonial)
# summary_2 <- get_bounds_unrest(obs)
# stats <- list(n = summary_1$n,
#               b_OLS = summary_1$b_OLS,
#               se_OLS = summary_1$se_OLS,
#               b_IV = summary_1$b_IV,
#               se_IV = summary_1$se_IV,
#               k_lower = summary_2$k$Lower,
#               r_uz_bound = ifelse(abs(summary_2$r_uz$Lower) == 1,
#                                   summary_2$r_uz$Upper,
#                                   summary_2$r_uz$Lower))
# make_I(stats, "Colonial Origins")
#
# # Part II of table
# ## First prior
# bounds_1 <- draw_bounds("logpgp95", "avexpr", "logem4", colonial, controls = NULL,
#                       r_TstarU_restriction = c(0, 0.9),
#                       k_restriction = c(0.00001, 0.6),
#                       n_draws = 5000, Jeffreys = TRUE)
# freq_1 <- summarize_bounds(bounds_1)
#
# posterior_1 <- draw_posterior("logpgp95", "avexpr", "logem4", colonial, controls = NULL,
#                              r_TstarU_restriction = c(0, 0.9),
#                             k_restriction = c(0.000001, 0.6),
#                             n_RF_draws = 1000,
#                           n_IS_draws = 1000,
#                           Jeffreys = TRUE,
#                           resample = FALSE)
# bayes_1 <- summarize_posterior(posterior_1)
#
# stats <- list(p_empty = freq_1$p_empty,
#               p_valid = freq_1$p_valid,
#               beta_lower_median = freq_1$restricted$median,
#               beta_upper_median = freq_1$restricted$median,
#               beta_lower_lower_bound = freq_1$restricted$lower[1],
#               beta_lower_upper_bound = freq_1$restricted$upper[1],
#               beta_upper_lower_bound = freq_1$restricted$lower[2],
#               beta_upper_upper_bound = freq_1$restricted$upper[2],
#               r_uz_median = bayes_1$HPDI$median[1],
#               beta_bayes_median = bayes_1$HPDI$median[2],
#               r_uz_lower_bound = bayes_1$HPDI$lower[1],
#               beta_bayes_lower_bound = bayes_1$HPDI$lower[2],
#               r_uz_upper_bound = bayes_1$HPDI$upper[1],
#               beta_bayes_upper_bound = bayes_1$HPDI$upper[2])
# make_II_III(stats, "$(\\kappa, \\rho_{T^*u}) \\in (0, 0.6] \\times [0, 0.9]")
#
# ## Second prior
# bounds_2 <- draw_bounds("logpgp95", "avexpr", "logem4", colonial, controls = NULL,
#                         r_TstarU_restriction = c(0, 0.9),
#                         k_restriction = c(0.6, 1),
#                         n_draws = 5000, Jeffreys = TRUE)
# freq_2 <- summarize_bounds(bounds_2)
#
# posterior_2 <- draw_posterior("logpgp95", "avexpr", "logem4", colonial, controls = NULL,
#                               r_TstarU_restriction = c(0, 0.9),
#                               k_restriction = c(0.6, 1),
#                               n_RF_draws = 1000,
#                               n_IS_draws = 1000,
#                               Jeffreys = TRUE,
#                               resample = FALSE)
# bayes_2 <- summarize_posterior(posterior_2)
#
# stats <- list(p_empty = freq_2$p_empty,
#               p_valid = freq_2$p_valid,
#               beta_lower_median = freq_2$restricted$median,
#               beta_upper_median = freq_2$restricted$median,
#               beta_lower_lower_bound = freq_2$restricted$lower[1],
#               beta_lower_upper_bound = freq_2$restricted$upper[1],
#               beta_upper_lower_bound = freq_2$restricted$lower[2],
#               beta_upper_upper_bound = freq_2$restricted$upper[2],
#               r_uz_median = bayes_2$HPDI$median[1],
#               beta_bayes_median = bayes_2$HPDI$median[2],
#               r_uz_lower_bound = bayes_2$HPDI$lower[1],
#               beta_bayes_lower_bound = bayes_2$HPDI$lower[2],
#               r_uz_upper_bound = bayes_2$HPDI$upper[1],
#               beta_bayes_upper_bound = bayes_2$HPDI$upper[2])
# make_II_III(stats, "$(\\kappa, \\rho_{T^*u}) \\in (0.6, 1] \\times [0, 0.9]")
#
#
# y_name <- "logpgp95"
# T_name <- "avexpr"
# z_name <- "logem4"
# data <- colonial
# controls <- NULL
# r_TstarU_restriction <- c(0, 0.9)
# k_restriction <- c(0.00001, 0.6)
# n_RF_draws <- 1000
# n_IS_draws <- 1000
# Jeffreys <- TRUE
# resample <- FALSE
#
# # r_TstarU_lower <- r_TstarU_min[!empty]
# # r_TstarU_upper <- r_TstarU_max[!empty]
# # k_lower <- k_min[!empty]
# # k_upper <- k_max[!empty]
# # obs <- obs_draws[!empty, ]
#
# ## Was Weber Wrong?
# load("./data/weber.rda")
# controls <- c("f_young", "f_jew", "f_fem", "f_ortsgeb", "f_pruss", "hhsize",
#               "lnpop", "gpop", "f_miss")
# summary_1 <- get_estimates("f_rw", "f_prot", "kmwittenberg", controls = controls, weber)
# obs <- get_observables("f_rw", "f_prot", "kmwittenberg", weber, controls = controls)
# summary_2 <- get_bounds_unrest(obs)
# stats <- list(n = summary_1$n,
#               b_OLS = summary_1$b_OLS,
#               se_OLS = summary_1$se_OLS,
#               b_IV = summary_1$b_IV,
#               se_IV = summary_1$se_IV,
#               k_lower = summary_2$k$Lower,
#               r_uz_bound = ifelse(abs(summary_2$r_uz$Lower) == 1,
#                                   summary_2$r_uz$Upper,
#                                   summary_2$r_uz$Lower))
# # make_I(stats, "Was Weber Wrong?")
