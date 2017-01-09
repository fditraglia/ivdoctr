# # Generating the tables for Colonial Origins
# load("./data/colonial.rda")
# ## Summary Stats
# get_estimates("logpgp95", "avexpr", "logem4", colonial)
# obs <- get_observables("logpgp95", "avexpr", "logem4", colonial)
# get_k_bounds_unrest(obs, tilde = FALSE)
# get_r_uz_bounds_unrest(obs)
# 
# draw_bounds("logpgp95", "avexpr", "logem4", colonial, controls = NULL, 
#             r_TstarU_restriction = c(0, 0.9), k_restriction = c(0, 0.6))
# draw_bounds("logpgp95", "avexpr", "logem4", colonial, controls = NULL, 
#             r_TstarU_restriction = c(0, 0.9), k_restriction = c(0, 0.6), 
#             Jeffreys = TRUE)
# 
# y_name <- "logpgp95"
# T_name <- "avexpr"
# z_name <- "logem4"
# data <- colonial
# controls <- NULL
# r_TstarU_restriction <- c(0, 0.9)
# k_restriction <- c(0, 0.6)
# n_draws <- 5000
# Jeffreys <- TRUE
# 
# # Generating the table for Weber
# load("./data/weber.rda")
# controls <- c("f_young", "f_jew", "f_fem", "f_ortsgeb", "f_pruss", "hhsize", 
#               "lnpop", "gpop", "f_miss")
# ## Summary stats
# get_estimates("f_rw", "f_prot", "kmwittenberg", controls = controls, weber)
# obs <- get_observables("f_rw", "f_prot", "kmwittenberg", weber, controls = controls)
# get_k_bounds_unrest(obs, tilde = TRUE)
# get_k_bounds_unrest(obs, tilde = FALSE)
# get_r_uz_bounds_unrest(obs)
