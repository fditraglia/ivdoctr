#' Rounds x to two decimal places
#' @param x Number to be rounded
#' @return Number rounded to 2 decimal places
myformat <- function(x){
  x <- ifelse(is.na(x), NA, format(round(x, 2), n_digits = 2, nsmall = 2))
}

#' Creates LaTeX code for parameter estimates
#' @param est Number
#' @return LaTeX string for the number
format_est <- function(est) {
  paste0('$', myformat(est), '$')
}

#' Creates LaTeX code for the standard error
#' @param se Standard error
#' @return LaTeX string for the standard error
format_se <- function(se) {
  paste0('$(', myformat(se), ')$')
}

#' Creates LaTeX code for the HPDI
#' @param bounds 2-element vector of the upper and lower HPDI bounds
#' @return LaTeX string of the HPDI
format_HPDI <- function(bounds) {
  paste0('$[', myformat(bounds[1]), ',', myformat(bounds[2]), ']$')
}

#' Makes LaTeX code to make a row of a table and shift by some amount of columns
#'   if necessary
#' @param char_vec Vector of characters to be collapsed into a LaTeX table
#' @param shift Number of columns to shift over
#' @return LaTeX string of the whole row of the table
make_tex_row <- function(char_vec, shift = 0) {
  out <- paste0(char_vec, collapse = ' & ')
  if (identical(shift, 0)) {
    out <- paste(out, '\\\\')
  } else {
    out <- paste(paste0(rep('&', shift), collapse = ''), out, '\\\\')
  }
  return(out)
}

#' Takes the OLS and IV estimates and converts it to a row of the LaTeX table
#' @param stats List with OLS and IV estimates and the bounds on kappa and r_uz
#' @param example_name Character string detailing the example
#' @return LaTeX code for the reduced-form estimates part of the table
make_I <- function(stats, example_name) {
  example_name <- paste(example_name, paste0('($n=', stats$n, '$)'))
  est <- with(stats, sapply(c(b_OLS, b_IV, k_lower, r_uz_bound), format_est))
  est <- make_tex_row(c(example_name, est))
  se <- with(stats, sapply(c(se_OLS, se_IV), format_se))
  se <- make_tex_row(se, shift = 1)
  paste(est, se, sep = '\\ ')
}

#' Takes frequentist or Bayesian simulations of beta and instrument validity
#'   and formats it into a LaTeX table
#' @param stats List containing parameter estimates and bounds to be output into table
#' @param prior_name String containing information on the user beliefs over
#'     endogeneity and measurement error
#' @return LaTeX string outputting a row of a table for those user restrictions
make_II_III <- function(stats, prior_name) {
  probs <- with(stats, sapply(c(p_empty, p_valid), format_est))

  medians <- with(stats, sapply(c(beta_lower,
                                  beta_upper,
                                  r_uz_median,
                                  beta_bayes_median), format_est))

  HPDIs <- with(stats, apply(matrix(c(r_uz_lower_bound, r_uz_upper_bound,
                                      beta_bayes_lower_bound, beta_bayes_upper_bound),
                                    ncol = 2, byrow = TRUE), 1, format_HPDI))

  row1 <- paste('\\hspace{2em}', prior_name, make_tex_row(c(probs, medians), shift = 5))
  row2 <- make_tex_row(HPDIs, shift = 9)
  paste(row1, row2, sep = '\\ ')
}

#' Generates table of parameter estimates given user restrictions and data
#'
#' @param y_name Character string with the column name of the dependent variable
#' @param T_name Character string with the column name of the endogenous regressor(s)
#' @param z_name Character string with the column name of the instrument(s)
#' @param data Data frame
#' @param controls Vector of character strings specifying the exogenous variables
#' @param robust Indicator for heteroskedasticity-robust standard errors
#' @param r_TstarU_restriction Matrix of desired mins and maxes for r_TstarU (must be same dimensions as k_restriction)
#' @param k_restriction Matrix of desired minx and maxes for kappa (must be same dimensions as r_TstarU_restriction)
#' @param n_draws Number of draws when generating frequentist-friendly draws of the covariance matrix
#' @param n_RF_draws Number of reduced-form draws
#' @param n_IS_draws Number of draws on the identified set
#' @param resample Indicator of whether or not to resample using magnification factor
#' @param Jeffreys Indicator of whether to draw covariance matrix using the Jeffreys prior to ensure positive definiteness
#' @param example_name Character string describing the example
#' @export
makeExample <- function(y_name, T_name, z_name, data, controls = NULL,
                        robust = FALSE, r_TstarU_restriction = NULL,
                        k_restriction = NULL, n_draws = 5000, n_RF_draws = 1000,
                        n_IS_draws = 1000, resample = FALSE, Jeffreys = FALSE,
                        example_name) {
  if (is.null(r_TstarU_restriction)) {
    r_TstarU_restriction <- matrix(c(-1, 1), nrow = 1)
  }
  if (is.null(k_restriction)) {
    k_restriction <- matrix(c(0, 1), nrow = 1)
  }
  if (nrow(r_TstarU_restriction) != nrow(k_restriction)) {
    if (is.null(r_TstarU_restriction)) {
      r_TstarU_restriction <- matrix(c(-1, 1), nrow = nrow(k_restriction), byrow = TRUE)
    } else if (is.null(k_restriction)) {
      k_restriction <- matrix(c(0, 1), nrow = nrow(r_TstarU_restriction), byrow = TRUE)
    } else {
      stop("Dimension mismatch between r_TstarU_restriction and k_restriction.
           Please make sure that if there are restrictions on kappa or r_TstarU
           that they are of the same dimension so that all examples are accounted
           for.")
    }
    }
  summary_stats <- get_estimates(y_name, T_name, z_name, data, controls, robust)
  obs <- get_observables(y_name, T_name, z_name, data, controls)
  bounds_unrest <- get_bounds_unrest(obs)

  stats_I <- list(n = summary_stats$n,
                  b_OLS = summary_stats$b_OLS,
                  se_OLS = summary_stats$se_OLS,
                  b_IV = summary_stats$b_IV,
                  se_IV = summary_stats$se_IV,
                  k_lower = bounds_unrest$k$Lower,
                  r_uz_bound = ifelse(abs(bounds_unrest$r_uz$Lower) == 1,
                                      bounds_unrest$r_uz$Upper,
                                      bounds_unrest$r_uz$Lower))
  headline <- make_I(stats_I, example_name)
  nExamples <- nrow(r_TstarU_restriction)
  exampleTex <- NULL
  for (i in 1:nExamples) {
    bounds <- draw_bounds(y_name, T_name, z_name, data, controls,
                          r_TstarU_restriction[i, ], k_restriction[i, ],
                          n_draws, Jeffreys)
    freq <- summarize_bounds(bounds)
    posterior <- draw_posterior(y_name, T_name, z_name, data, controls,
                                r_TstarU_restriction[i, ], k_restriction[i, ],
                                n_RF_draws, n_IS_draws, Jeffreys, resample)
    bayes <- summarize_posterior(posterior)

    # Compute covering beta interval
    center <- bounds$center
    beta_bounds <- cbind(bounds$restricted$beta_lower, bounds$restricted$beta_upper)
    beta_interval <- getInterval(beta_bounds, center)
    stats <- list(p_empty = freq$p_empty,
                  p_valid = freq$p_valid,
                  beta_lower = min(beta_interval),
                  beta_upper = max(beta_interval),
                  r_uz_median = bayes$HPDI$median[1],
                  beta_bayes_median = bayes$HPDI$median[2],
                  r_uz_lower_bound = bayes$HPDI$lower[1],
                  beta_bayes_lower_bound = bayes$HPDI$lower[2],
                  r_uz_upper_bound = bayes$HPDI$upper[1],
                  beta_bayes_upper_bound = bayes$HPDI$upper[2])
    newRow <- make_II_III(stats, paste0("$(\\kappa, \\rho_{T^*u}) \\in (",
                                        ifelse(k_restriction[i, 1] < 0.01,
                                               0, k_restriction[i, 1]), ",",
                                        k_restriction[i, 2],
                                        "] \\times [",
                                        r_TstarU_restriction[i, 1], ",",
                                        r_TstarU_restriction[i, 2], "]$"))
    exampleTex <- paste(exampleTex, newRow, sep = "\\ ")
  }
  output <- paste(headline, exampleTex, sep = "\\ ")
  return(output)
    }

#' Generates header LaTeX code for continuous table
#'
table_header_fn <- function() {
  "\\begin{tabular}{lcccccccccc}
  \\hline
  \\hline
  &\\multicolumn{4}{c}{(I) Summary Statistics}
  &\\multicolumn{4}{c}{(II) Frequentist-Friendly}
  &\\multicolumn{2}{c}{(III) Full Bayesian} \\\\
  \\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-11}
  & OLS & IV & $\\underline{\\kappa}$ & $\\underline{\\tilde{\\rho}}_{uz}/\\tilde{\\bar{\\rho}}_{uz}$ & $\\mathbb{P}(\\varnothing)$ & $\\mathbb{P}(\\mbox{Valid})$ & $\\underline{\\beta}$ & $\\bar{\\beta}$ & $\\rho_{uz}$ & $\\beta$ \\\\
  \\\\"
}

#' Generates footer LaTeX code for continuous table
#'
table_footer_fn <- function() {
  "\\hline
  \\end{tabular}"
}

#' Generates table of parameter estimates given user restrictions and data
#'
#' @param file Character string with the path name of the .tex file to be saved
#' @param ... Arguments of TeX code for individual examples to be combined
#'   into a single table
#' @export
makeTable <- function(file, ...) {

  table_examples <- c()
  list_examples <- list(...)

  for (i in 1:length(list(...))) {
    table_examples <- c(table_examples, list_examples[[i]], "\\\\")
  }

  cat(table_header_fn(), '\\\\', table_examples, table_footer_fn(),
      sep = '\n', file = file)
}
