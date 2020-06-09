#' @import grDevices
#' @import graphics
#' @import rgl
#'
NULL

# Generates plot for a given restriction on r_TstarU and kappa
#' Plot ivdoctr Restrictions
#'
#' @param r_TstarU_restriction Matrix of desired mins and maxes for r_TstarU (must be same dimensions as k_restriction)
#' @param k_restriction Matrix of desired minx and maxes for kappa (must be same dimensions as r_TstarU_restriction)
#' @param n_grid Number of points to put in grid
#' @param n_colors Number of colors to use
#' @param y_name Character string with the column name of the dependent variable
#' @param T_name Character string with the column name of the endogenous regressor(s)
#' @param z_name Character string with the column name of the instrument(s)
#' @param data Data frame
#' @param controls Vector of character strings specifying the exogenous variables
#' @param gray_k Matrix of kappa restrictions to recolor graph as gray
#' @param gray_rTstarU Matrix of rTstarU restrictions to recolor graph as gray
#' @param fence Vector of left, bottom, right, and top corners of rectangle
#' @param theta Graphing parameters for orienting plot
#' @param phi Graphing parameters for orienting plot
#'
#' @return Interactive 3d plot which can be oriented and saved using rgl.snapshot()
#' @export
#'
plot_3d_beta <- function(y_name, T_name, z_name, data, controls = NULL,
                         r_TstarU_restriction = NULL, k_restriction = NULL,
                         n_grid = 30, n_colors = 500, fence = NULL,
                         gray_k = NULL, gray_rTstarU = NULL,
                         theta = 0, phi = 15) {
  color <- NULL
  obs <- get_observables(y_name, T_name, z_name, data, controls)

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
    r_TstarU_lower <- min(r_TstarU_restriction[i, ])
    r_TstarU_upper <- max(r_TstarU_restriction[i, ])
    r_TstarU <- seq(r_TstarU_lower, r_TstarU_upper, length.out = n_grid)

    # We will plot k, but get_r_uz and get_beta both take k_tilde as inputs
    k <- seq(k_lower[i], k_upper[i], length.out = n_grid)
    k_tilde <- (k - obs$T_Rsq) / (1 - obs$T_Rsq)
    fullGrid <- expand.grid(r_TstarU = r_TstarU, k_tilde = k_tilde)
    r_uz <- get_r_uz(fullGrid$r_TstarU, fullGrid$k_tilde, obs)
    beta <- get_beta(fullGrid$r_TstarU, fullGrid$k_tilde, obs)

    # Creating grid of points and removing infinite values of beta
    graphData <- data.table(cbind(fullGrid, r_uz = r_uz, beta = beta))[beta < Inf & beta > -Inf]

    # Generating colors for beta. Red is negative, blue is positive, and gray is
    # based on a separate argument
    graphData[beta >= 0, "color" := "blue"]
    graphData[beta < 0, "color" := "red"]

    blues <- colorRampPalette(c("white", "blue"))(n_colors / 2)
    reds <- colorRampPalette(c("red", "white"))(n_colors / 2)
    gray <- colorRampPalette("gray")(1)

    if (!is.null(gray_k)) {
      graphData[k_tilde < min(gray_k[i, ]), ':=' (color = "gray", hexCol = gray)]
      graphData[k_tilde > max(gray_k[i, ]), ':=' (color = "gray", hexCol = gray)]
    }
    if (!is.null(gray_rTstarU)) {
      graphData[r_TstarU < min(gray_rTstarU[i, ]), ':=' (color = "gray", hexCol = gray)]
      graphData[r_TstarU > max(gray_rTstarU[i, ]), ':=' (color = "gray", hexCol = gray)]
    }

    blueCol <- map2color(graphData[color == "blue"]$beta, blues)
    redCol <- map2color(graphData[color == "red"]$beta, reds)

    graphData[color == "blue", "hexCol" := blueCol]
    graphData[color == "red", "hexCol" := redCol]

    # Rearranging data to match expected arguments for persp3d()
    rho_tu <- matrix(graphData$r_TstarU, ncol = n_grid)
    kappa <- matrix(graphData$k_tilde, ncol = n_grid)
    rho_uz <- matrix(graphData$r_uz, ncol = n_grid)

    view3d(theta, phi)
    persp3d(x = rho_tu, y = kappa, z = rho_uz, col = graphData$hexCol,
            xlim = c(r_TstarU_lower, r_TstarU_upper), ylim = c(k_lower, k_upper),
            xlab = "Rho_tu", ylab = "Kappa", zlab = "Rho_uz")

    if (!is.null(fence)) {
      left <- r_TstarU[which.min(abs(r_TstarU - fence[1]))]
      bottom <- k[which.min(abs(k - fence[2]))]
      right <- r_TstarU[which.min(abs(r_TstarU - fence[3]))]
      top <- k[which.min(abs(k - fence[4]))]
      k_diff <- diff(k)[1]
      r_TstarU_diff <- diff(r_TstarU)[1]

      my_rect <- rect_points(left, bottom, right, top, r_TstarU_diff, k_diff)
      my_rect$z <- get_r_uz(my_rect$x, my_rect$y, obs)
      lines3d(my_rect)
    }
  }
}

#' Generates a custom color palette given a vector of numbers
#'
#' @param x Vector of numbers
#' @param pal Palette function generate from colorRampPalette
#' @param limits Limits on the numeric sequence
#'
#' @return Hex values for colors
#'
map2color <- function(x, pal, limits = NULL){
  if(is.null(limits)) limits = range(x)
  pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1), all.inside = TRUE)]
}
