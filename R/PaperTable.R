# # ---------------------------- Table header and footer -------------------------
# table_header <- "\\begin{tabular}{lcccccccccc}
#   \\hline
#   \\hline
#   &\\multicolumn{4}{c}{(I) Summary Statistics}
#   &\\multicolumn{4}{c}{(II) Frequentist-Friendly}
#   &\\multicolumn{2}{c}{(III) Full Bayesian} \\\\
#   \\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-11}
#   & OLS & IV & $\\underline{\\kappa}$ & $\\underline{\\rho}_{uz}/\\bar{\\rho}_{uz}$ & $\\mathbb{P}(\\varnothing)$ & $\\mathbb{P}(\\mbox{Valid})$ & $\\underline{\\beta}$ & $\\bar{\\beta}$ & $\\rho_{uz}$ & $\\beta$ \\\\
#   \\\\"
#
# table_footer <- "\\hline
#                  \\end{tabular}"
#
# # ----------------------------- Colonial Origins -------------------------------
# r_TstarU_restriction <- matrix(c(0, 0.9), nrow = 2, ncol = 2, byrow = TRUE)
# k_restriction <- matrix(c(0.0001, 0.6, 0.6, 1), nrow = 2, ncol = 2, byrow = TRUE)
#
# colonial_example <- makeExample("logpgp95", "avexpr", "logem4", colonial,
#                                 controls = NULL, robust = FALSE,
#                                 r_TstarU_restriction = r_TstarU_restriction,
#                                 k_restriction = k_restriction,
#                                 n_draws = 5000, n_RF_draws = 1000,
#                                 n_IS_draws = 1000, resample = TRUE,
#                                 Jeffreys = TRUE,
#                                 example_name = "Colonial Origins")
#
# # ---------------------------- Weber Example -----------------------------------
# r_TstarU_restriction <- matrix(c(-0.9, 0), nrow = 2, ncol = 2, byrow = TRUE)
# k_restriction <- matrix(c(0.0001, 1, 0.8, 1), nrow = 2, ncol = 2, byrow = TRUE)
#
# weber_example <- makeExample("f_rw", "f_prot", "kmwittenberg", weber,
#                                 controls = c("f_young", "f_jew", "f_fem",
#                                              "f_ortsgeb", "f_pruss", "hhsize",
#                                              "lnpop", "gpop", "f_miss",
#                                              "f_blind", "f_deaf", "f_dumb"),
#                              robust = FALSE,
#                              r_TstarU_restriction = r_TstarU_restriction,
#                              k_restriction = k_restriction, n_draws = 5000,
#                              n_RF_draws = 1000, n_IS_draws = 1000,
#                              resample = TRUE, Jeffreys = TRUE,
#                              example_name = "Was Weber Wrong?")
#
# cat(table_header, '\\\\',
#     colonial_example, '\\\\',
#     weber_example, '\\\\',
#     table_footer,
#     sep = '\n', file = "./data-raw/continuous_table.tex")
#
