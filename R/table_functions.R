# #### This file contains three functions:
# 1. texLine_header(): Initial tex line for each example
# 2. texLine_prior(): Line for each prior with each example
# 3. makeTexTable(): Carries out numerical calculations and uses the above
#    functions as building blocks.

# Note: Only makeTexTable() is directly used by other programs. Below are
# instructions on what type of input is required.

# ------------------------ Generic Helper Functions
myformat <- function(x){
  x <- format(round(x, 2), n_digits = 2, nsmall = 2)
}
format_est <- function(est) {
  paste0('$', myformat(est), '$')
}
format_se <- function(se) {
  paste0('$(', myformat(se), ')$')
}
format_HPDI <- function(lower, upper) {
  paste0('$[', myformat(lower), ',', myformat(upper), ']$')
}

make_tex_row <- function(char_vec, shift = 0) {
  out <- paste0(char_vec, collapse = ' & ')
  if (identical(shift, 0)) {
    out <- paste(out, '\\\\')
  } else {
    out <- paste(paste0(rep('&', shift), collapse = ''), out, '\\\\')
  }
  return(out)
}

# ------------------------ Functions Specific to Continuous Case
make_I <- function(stats, example_name) {
  example_name <- paste(example_name, paste0('($n=', stats$n, '$)'))
  est <- with(stats, sapply(c(b_OLS, b_IV, k_lower, r_uz_bound), format_est))
  est <- make_tex_row(c(example_name, est))
  se <- with(stats, sapply(c(se_OLS, se_IV), format_se))
  se <- make_tex_row(se, shift = 1)
  paste(est, se, sep = '\n')
}

make_II <- function(stats, prior_name) {
    probs <- with(stats, sapply(c(p_empty, p_valid), format_est))

    medians <- with(stats, sapply(c(beta_lower_median, beta_upper_median), format_est))

    HPDIs <- with(stats, sapply(c(beta_lower_lower_bound,
                                  beta_lower_upper_bound,
                                  beta_upper_lower_bound,
                                  beta_upper_upper_bound),
                                format_est))

  row1 <- paste('\\hspace{2em}', prior_name, make_tex_row(c(probs, medians), shift = 5))
  row2 <- make_tex_row(HPDIs, shift = 7)
  paste(row1, row2, sep = '\n')
}

make_III <- function(stats, prior_name) {
    medians <- c(format_est(stats['b_lower', 'median']),
                 format_est(stats['b_upper', 'median']))

    HPDIs <- c(format_HPDI(stats['b_lower', 'lower'],
                           stats['b_lower', 'upper']),
               format_HPDI(stats['b_upper', 'lower'],
                           stats['b_upper', 'upper']))

  row1 <- paste('\\hspace{2em}', prior_name, make_tex_row(medians, shift = 5))
  row2 <- make_tex_row(HPDIs, shift = 7)
  paste(row1, row2, sep = '\n')
}

# ------------------------ Actual Table Functions
## Creates the line in each example with the summary statistics
texLine_header <- function(title, n, b_ols, se_ols, b_iv, se_iv,
                           underline_kappa, RzuBound) {

  paste0(title, format_se(n), " & ", format(round(b_ols, 2), nsmall = 2),
         " & ", format(round(b_iv, 2), nsmall = 2), " & ", round(underline_kappa, 2), " & ", round(RzuBound, 2), " \\\\ & (",
        round(se_ols, 2), ") & (", round(se_iv, 2),") \\\\")
}

## Creates the line for one example and one prior for that example
texLine_prior <- function(textprior, prob_empty, prob_valid,
                          median_underline_beta, median_overline_beta,
                          median_rho_uz, median_beta, lower90_rho_uz, upper90_rho_uz,
                          lower90_beta, upper90_beta,
                          lower90_underline_beta, upper90_underline_beta,
                          lower90_overline_beta, upper90_overline_beta) {

  if (prob_empty > 0) {
    prob_empty <- round(prob_empty, 2)
    prob_valid <- "---"
    median_underline_beta <- "---"
    median_overline_beta <- "---"
    coverage90_underline_beta <- "\\mbox{---},\\mbox{---}"
    coverage90_overline_beta <- "\\mbox{---},\\mbox{---}"
    median_rho_uz <- "---"
    median_beta <- "---"
    coverage90_rho_uz <- "\\mbox{---},\\mbox{---}"
    coverage90_beta <- "\\mbox{---},\\mbox{---}"
  } else {
    prob_empty <- format(round(prob_empty, 2), nsmall = 2)
    prob_valid <- format(round(prob_valid, 2), nsmall = 2)
    median_underline_beta <- format(round(median_underline_beta, 2), nsmall = 2)
    median_overline_beta <- format(round(median_overline_beta, 2), nsmall = 2)
    coverage90_underline_beta <- paste(format(round(lower90_underline_beta, 2), nsmall = 2),",",format(round(upper90_underline_beta, 2), nsmall = 2))
    coverage90_overline_beta <- paste(format(round(lower90_overline_beta, 2), nsmall = 2),",",format(round(upper90_overline_beta, 2), nsmall = 2))
    median_rho_uz <- format(round(median_rho_uz, 2), nsmall = 2)
    median_beta <- format(round(median_beta, 2), nsmall = 2)
    coverage90_rho_uz <- paste(format(round(lower90_rho_uz, 2), nsmall = 2),",",format(round(upper90_rho_uz, 2), nsmall = 2))
    coverage90_beta <- paste(format(round(lower90_beta, 2), nsmall = 2),",",format(round(upper90_beta, 2), nsmall = 2))
  }

  paste("\\hspace{2em}",textprior, "&&&&&", prob_empty,"&", prob_valid ,
        "&", median_underline_beta ,"&", median_overline_beta ,"&", median_rho_uz, "&", median_beta,"\\\\
        &&&&&&& $[", coverage90_underline_beta, "]$ & $[", coverage90_overline_beta,"]$ & $[", coverage90_rho_uz ,"]$ & $[",coverage90_beta,"]$ \\\\",sep = "")
}

### Creates a String of Tex Code with the whole table.
# The input must be of the format:
# list( listExample1, listExample2, .... )
# The list for each example must be written as
# listExample1 <- list("Name of Example (shown in table)",database, ivdoctr_object1,ivdoctr_object2, ...)
# ivdoctr_object is the output of using the "ivdoctr" package.

makeTexTable <- function(x) {

  L <- length(x)

  table <- c()

  table_header <- "\\begin{tabular}{lcccccccccc}
  \\hline \\hline
  &\\multicolumn{4}{c}{(I) Summary Statistics}
  &\\multicolumn{4}{c}{(II) Frequentist-Friendly}
  &\\multicolumn{2}{c}{(III) Full Bayesian} \\\\
  \\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-11}
  & OLS & IV & $\\underline{\\kappa}$ & $\\underline{\\rho}_{uz}/\\bar{\\rho}_{uz}$ & $\\mathbb{P}(\\varnothing)$ & $\\mathbb{P}(\\mbox{Valid})$ & $\\underline{\\beta}$ & $\\bar{\\beta}$ & $\\rho_{uz}$ & $\\beta$ \\\\
  \\\\"

  table <- cbind(table,table_header)

  for (i in 1:L) {

    example_input <- x[[i]]
    example_output <- c()

    title <- example_input[[1]]
    data <- example_input[[2]]
    seednumber <- example_input[[3]]
    ivdoctr_output <- eval(example_input[[4]])

    ## print(title)
    ## print(data)

    n <- ivdoctr_output$summary$n
    b_ols <- coef(summary(lm(reformulate(c(ivdoctr_output$summary$x_name,ivdoctr_output$summary$controls), response = ivdoctr_output$summary$y_name), data)))[ivdoctr_output$summary$x_name,"Estimate"]
    se_ols <- coef(summary(lm(reformulate(c(ivdoctr_output$summary$x_name,ivdoctr_output$summary$controls), response = ivdoctr_output$summary$y_name), data)))[ivdoctr_output$summary$x_name,"Std. Error"]
    b_iv <- coef(summary(ivreg(formula=reformulate(c(ivdoctr_output$summary$x_name,ivdoctr_output$summary$controls), response = ivdoctr_output$summary$y_name), instruments=reformulate(c(ivdoctr_output$summary$z_name,ivdoctr_output$summary$controls)),data=data)))[ivdoctr_output$summary$x_name,"Estimate"]
    se_iv <- coef(summary(ivreg(formula=reformulate(c(ivdoctr_output$summary$x_name,ivdoctr_output$summary$controls), response = ivdoctr_output$summary$y_name), instruments=reformulate(c(ivdoctr_output$summary$z_name,ivdoctr_output$summary$controls)),data=data)))[ivdoctr_output$summary$x_name,"Std. Error"]

    underline_kappa <- toKappa(getUnderlineKappa(ivdoctr_output$MLE$Sigma),ivdoctr_output$MLE$xRsq)
    RzuBound <- ifelse(getRzuBounds_noprior(ivdoctr_output$MLE$Sigma,underline_kappa)[2]==1,getRzuBounds_noprior(ivdoctr_output$MLE$Sigma,underline_kappa)[1],getRzuBounds_noprior(ivdoctr_output$MLE$Sigma,underline_kappa)[2])

    example_output <- cbind(example_output,texLine_header(title,n,b_ols,se_ols,b_iv,se_iv,underline_kappa,RzuBound))
    rm(ivdoctr_output)

    for(j in 4:length(example_input)) {

      print(paste("Processing ...",title, "... Prior ...",j-3," / ", length(example_input)-3 ))

      set.seed(seednumber)
      ivdoctr_output <- eval(example_input[[j]])
      textprior <- paste("$(\\kappa, \\rho_{T^*u}) \\in (",
                         ivdoctr_output$summary$prior$K[1], ",", ivdoctr_output$summary$prior$K[2], "]\\times [",
                         ivdoctr_output$summary$prior$Rxsu[1], ",", ivdoctr_output$summary$prior$Rxsu[2], "]$")

      prob_empty <- 1-mean(ivdoctr_output$sim$diagnostic$Klower < ivdoctr_output$summary$prior$K[2])

      if( prob_empty == 0 ) {

        prob_valid <- sum((sapply(ivdoctr_output$sim$draws_list$Rzu, function(x) all(!is.na(x)))) &
                            (sapply(ivdoctr_output$sim$draws_list$Rzu, min) < 0) &
                            (sapply(ivdoctr_output$sim$draws_list$Rzu, max) > 0)) /
          length(ivdoctr_output$sim$Sigma)

        beta_underline <- sapply(ivdoctr_output$sim$Sigma,function(x) getBeta_min(x,ivdoctr_output$summary$prior$Rxsu,ivdoctr_output$summary$prior$K))
        median_underline_beta <- median(beta_underline)
        coverageinterval_underline_beta <- HPDinterval(as.mcmc(beta_underline), prob = 0.9)
        lower90_underline_beta <- coverageinterval_underline_beta[1]
        upper90_underline_beta <- coverageinterval_underline_beta[2]

        beta_overline <- sapply(ivdoctr_output$sim$Sigma,function(x) getBeta_max(x,ivdoctr_output$summary$prior$Rxsu,ivdoctr_output$summary$prior$K))
        median_overline_beta <- median(beta_overline)
        coverageinterval_overline_beta <- HPDinterval(as.mcmc(beta_overline), prob = 0.9)
        lower90_overline_beta <- coverageinterval_overline_beta[1]
        upper90_overline_beta <- coverageinterval_overline_beta[2]

        median_rho_uz <- median(ivdoctr_output$sim$draws$Rzu)
        median_beta <- median(ivdoctr_output$sim$draws$Beta)

        coverageinterval_rho_uz <- HPDinterval(as.mcmc(ivdoctr_output$sim$draws$Rzu), prob = 0.9)
        lower90_rho_uz <- coverageinterval_rho_uz[1]
        upper90_rho_uz <- coverageinterval_rho_uz[2]

        coverageinterval_beta <- HPDinterval(as.mcmc(ivdoctr_output$sim$draws$Beta), prob = 0.9)
        lower90_beta <- coverageinterval_beta[1]
        upper90_beta <- coverageinterval_beta[2]

        example_output <- cbind(example_output,texLine_prior(textprior,prob_empty,
                                                             prob_valid,median_underline_beta, median_overline_beta,
                                                             median_rho_uz,median_beta,lower90_rho_uz,upper90_rho_uz,
                                                             lower90_beta,upper90_beta,
                                                             lower90_underline_beta,upper90_underline_beta,
                                                             lower90_overline_beta,upper90_overline_beta))

      } else {
        example_output <- cbind(example_output,texLine_prior(textprior,prob_empty,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0,0,0,0))

      }

      rm(ivdoctr_output)

    }

    example_output <- cbind(example_output, "\\\\")
    table <- cbind(table,example_output)
  }

  table <- cbind(table,"\\hline
                 \\end{tabular}")

  return(table)

}


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
