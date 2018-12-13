#' @importFrom AER ivreg
#' @importFrom coda HPDinterval as.mcmc
#' @importFrom sandwich vcovHC
#' @import grDevices
#' @import graphics
NULL

#' Computes OLS and IV estimates
#' @param y_name Character vector of the name of the dependent variable
#' @param T_name Character vector of the names of the preferred regressors
#' @param z_name Character vector of the names of the instrumental variables
#' @param data Data to be analyzed
#' @param controls Character vector containing the names of the exogenous regressors
#' @param robust Boolean of whether to compute heterskedasticity-robust standard errors
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
  p_valid <- get_p_valid(draws)
  list(r_TstarU_restriction = draws$r_TstarU_restriction,
       k_restriction = draws$k_restriction,
       p_empty = mean(draws$empty),
       p_valid = p_valid,
       HPDI = HPDI)
}


# Generates plot for a given restriction on r_TstarU and kappa
plot_3d_beta <- function(obs, r_TstarU_restriction = NULL, k_restriction = NULL,
                         n_grid = 30, n_colors = 500, fence = NULL,
                         make_gray = function(r_TstarU, k) TRUE) {

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

  if (!is.null(k_restriction)) {
    k_lower <- pmax(get_bounds_unrest(obs)$k$Lower, min(k_restriction))
    k_upper <- pmin(1, max(k_restriction))
  }

  nPlots <- max(nrow(k_restriction), nrow(r_TstarU_restriction))
  plots <- list()
  for (i in nPlots) {
    r_TstarU_lower <- r_TstarU_restriction[i, 1]
    r_TstarU_upper <- r_TstarU_restriction[i, 2]
    r_TstarU <- seq(r_TstarU_lower, r_TstarU_upper, length.out = n_grid)
    # We will plot k, but get_r_uz and get_beta both take k_tilde as inputs
    k <- seq(k_lower[i], k_upper[i], length.out = n_grid)
    k_tilde <- (k - obs$T_Rsq) / (1 - obs$T_Rsq)
    fullGrid <- expand.grid(r_TstarU, k_tilde)
    r_uz <- get_r_uz(fullGrid[, 1], fullGrid[, 2], obs)
    beta <- get_beta(fullGrid[, 1], fullGrid[, 2], obs)

    is_gray <- outer(r_TstarU, k, make_gray)
    is_gray_facet <- 1 < (is_gray[-1, -1] + is_gray[-1, -n_grid] +
                          is_gray[-n_grid, -1] + is_gray[-n_grid, -n_grid])
    beta_facet <- 0.25 * (beta[-1, -1] + beta[-1, -n_grid] +
                          beta[-n_grid, -1] + beta[-n_grid, -n_grid])

    neg <- (beta_facet < 0) & (!is_gray_facet)
    pos <- (beta_facet >= 0) & (!is_gray_facet)
    n_neg <- ceiling(mean(neg) * n_colors)
    n_pos <- n_colors - n_neg
    neg_colors <- colorRampPalette(c('red', 'white'))(n_neg)
    pos_colors <- colorRampPalette(c('white', 'blue'))(n_pos)
    facet_colors <- matrix(NA, n_grid - 1, n_grid - 1)

    neg_color_bins <- cut(beta_facet[neg], n_neg)
    pos_color_bins <- cut(beta_facet[pos], n_pos)
    facet_colors[neg] <- neg_colors[neg_color_bins]
    facet_colors[pos] <- pos_colors[pos_color_bins]
    facet_colors[is_gray_facet] <- gray(0.6)

    bin_levels <- c(levels(neg_color_bins), levels(pos_color_bins))
    bins_lower <- as.numeric(sub("\\((.+),.*", "\\1", bin_levels))
    bins_upper <- as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", bin_levels))
    all_colors <- c(neg_colors, pos_colors)

    layout(matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(1, 4),
           heights = c(4, 4))

    plot(1, 1, t = 'n', ylim = c(bins_lower[1], bins_upper[n_colors]),
         ylab = '$\\beta$', xlim = c(0,1), xaxt = 'n', yaxt = 's', xlab = '',
         xaxs = 'i', yaxs = 'i')

    for (i in seq_len(n_colors)) {
        polygon(c(0, 0, 1, 1), c(bins_lower[i], bins_upper[i], bins_upper[i],
                              bins_lower[i]), col = all_colors[i], border = NA)
    }
    VT <- persp(r_TstarU, k, r_uz, ticktype = 'detailed', nticks = 4, phi = 30,
                theta = 225, xlab = '$\\rho_{T^*u}$', ylab = '$\\kappa$',
                zlab = '$\\rho_{zu}$', col = facet_colors, cex.axis = 0.75)

    if (!is.null(fence)) {

      left <- r_TstarU[which.min(abs(r_TstarU - fence[1]))]
      bottom <- k[which.min(abs(k - fence[2]))]
      right <- r_TstarU[which.min(abs(r_TstarU - fence[3]))]
      top <- k[which.min(abs(k - fence[4]))]
      k_diff <- diff(k)[1]
      r_TstarU_diff <- diff(r_TstarU)[1]

      my_rect <- rect_points(left, bottom, right, top, r_TstarU_diff, k_diff)
      my_rect$z <- get_r_uz(my_rect$x, my_rect$y, obs)
      my_rect <- with(my_rect, trans3d(x, y, z, VT))
      with(my_rect, points(x, y, type = 'l', lwd = 3))

    }
    par(mfrow = c(1, 1))
  }
}
