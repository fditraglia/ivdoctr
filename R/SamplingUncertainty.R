#' @import stats
#' @importFrom MASS mvrnorm
#' @importFrom utils head tail
NULL

#' Draws covariance matrix using the Jeffrey's Prior
#'
#' @param y Vector of dependent variable
#' @param Tobs Matrix containing data for the preferred regressor
#' @param z Matrix containing data for the instrumental variable
#' @param n_draws Integer number of draws to perform
#' @param k Number of covariates, including the intercept
#'
#' @return Array of covariance matrix draws
draw_sigma_jeffreys <- function(y, Tobs, z, k, n_draws) {
  n <- length(y)
  v <- n - k + 3 + 1
  S <- (n - 1) * cov(cbind(Tobs, y, z))
  Sigma_draws <- rinvwish(n_draws, v, S)
  rownames(Sigma_draws) <- colnames(Sigma_draws) <- c("Tobs", "y", "z")
  return(Sigma_draws)
}

#' Simulates different data draws
#'
#' This function takes the data and simulates potential draws of data from
#'   the properties of the observed data.
#'
#' @param y_name Character vector of the name of the dependent variable
#' @param T_name Character vector of the names of the preferred regressors
#' @param z_name Character vector of the names of the instrumental variables
#' @param data Data to be analyzed
#' @param controls Character vector containing the names of the exogenous regressors
#' @param n_draws Integer number of simulations to draw
#'
#' @return Data frame containing covariances, correlations, and R-squares for
#'   each data simulation
draw_observables <- function(y_name, T_name, z_name, data, controls = NULL,
                             n_draws = 5000) {

  # Project out control regressors if present
  if (!is.null(controls)) {
    y <- resid(lm(reformulate(controls, response = y_name), data))
    T_reg <- lm(reformulate(controls, response = T_name), data)
    z_reg <- lm(reformulate(controls, response = z_name), data)
    Tobs <- resid(T_reg)
    T_Rsq <- summary(T_reg)$r.squared
    z <- resid(z_reg)
    z_Rsq <- summary(z_reg)$r.squared
  } else {
    y <- get(y_name, data)
    Tobs <- get(T_name, data)
    z <- get(z_name, data)
    T_Rsq <- z_Rsq <- 0 # No controls is the same as controls that are
                        # uncorrelated with both Tobs and z
  }
  k <- length(controls) + 1
  Sigma_draws <- draw_sigma_jeffreys(y, Tobs, z, k, n_draws)
  not_positive_definite <- rep(FALSE, n_draws)

  s2_T <- Sigma_draws["Tobs", "Tobs", ]
  s2_y <- Sigma_draws["y", "y", ]
  s2_z <- Sigma_draws["z", "z", ]
  s_T <- sqrt(s2_T)
  s_y <- sqrt(s2_y)
  s_z <- sqrt(s2_z)
  s_Ty <- Sigma_draws["Tobs", "y", ]
  s_Tz <- Sigma_draws["Tobs", "z", ]
  s_zy <- Sigma_draws["z", "y", ]
  r_Ty <- s_Ty / (s_T * s_y)
  r_Tz <- s_Tz / (s_T * s_z)
  r_zy <- s_zy / (s_z * s_y)

  data.frame(n = rep(nrow(data), n_draws),
             T_Rsq = rep(T_Rsq, n_draws),
             z_Rsq = rep(z_Rsq, n_draws),
             s2_T = s2_T,
             s2_y = s2_y,
             s2_z = s2_z,
             s_Ty = s_Ty,
             s_Tz = s_Tz,
             s_zy = s_zy,
             r_Ty = r_Ty,
             r_Tz = r_Tz,
             r_zy = r_zy,
             not_positive_definite = not_positive_definite)
}

#' Computes bounds for simulated data
#'
#' This function takes data and user restrictions on measurement error and
#'   endogeneity and simulates data and the resulting bounds on instrument
#'   validity.
#' @param y_name Character vector of the name of the dependent variable
#' @param T_name Character vector of the names of the preferred regressors
#' @param z_name Character vector of the names of the instrumental variables
#' @param data Data to be analyzed
#' @param controls Character vector containing the names of the exogenous regressors
#' @param r_TstarU_restriction 2 element vector containing the min and max
#'   imposed on r_TstarU
#' @param k_restriction 2-element vector containing the min and max imposed on kappa
#' @param n_draws Integer number of simulations to draw
#'
#' @return List containing simulated data observables (covariances,
#'   correlations, and R-squares), indications of whether the identified set
#'   is empty, the unrestricted and restricted bounds on instrumental relevance,
#'   instrumental validity, and measurement error.
draw_bounds <- function(y_name, T_name, z_name, data, controls = NULL,
                        r_TstarU_restriction = NULL, k_restriction = NULL,
                        n_draws = 5000) {

  obs_draws <- draw_observables(y_name, T_name, z_name, data, controls, n_draws)
  n_draws <- n_draws + 1 # adding average draw at the end.
  obs_draws <- rbind(obs_draws, lapply(obs_draws, mean)) # making final row the mean of observables
  unrestricted_bounds <- get_bounds_unrest(obs_draws)
  empty <- NULL
  k_tilde_lower <- unrestricted_bounds$k_tilde$Lower

  if (!is.null(r_TstarU_restriction) | !is.null(k_restriction)) {
    if (!is.null(r_TstarU_restriction)) {
      r_TstarU_min <- rep(min(r_TstarU_restriction), n_draws)
      r_TstarU_max <- rep(max(r_TstarU_restriction), n_draws)
    }

    if (!is.null(k_restriction)) {
      # User states beliefs over kappa but we work with kappa_tilde
      k_min <- with(obs_draws, (min(k_restriction) - T_Rsq) / (1 - T_Rsq))
      k_max <- with(obs_draws, (max(k_restriction) - T_Rsq) / (1 - T_Rsq))
      k_min <- pmax(k_tilde_lower, k_min) # vector: could vary with obs_draws row
      k_max <- pmin(1, k_max)
    } else {
      k_min <- k_tilde_lower # vector: varies with obs_draws row
      k_max <- rep(1, n_draws)
    }

    # We ensure above that but the user may have specified a k_max that is less
    # than some elements of k_tilde_lower as in one of the examples for Colonial
    # Origins from the paper. When this occurs, the identified set is empty.
    empty <- k_max < k_tilde_lower

    # Only compute the bounds for the non-empty identified sets
    beta_lower <- rep(NA, n_draws)
    beta_lower[which(!empty)] <- get_beta_lower(r_TstarU_max[!empty], k_min[!empty],
                                                k_max[!empty], obs_draws[!empty, ])
    beta_upper <- rep(NA, n_draws)
    beta_upper[which(!empty)] <- get_beta_upper(r_TstarU_min[!empty], k_min[!empty],
                                                k_max[!empty], obs_draws[!empty, ])
    r_uz_restricted <- data.frame(min = rep(NA, n_draws),
                                  max = rep(NA, n_draws))
    r_uz_restricted[which(!empty), ] <- get_r_uz_bounds(r_TstarU_min[!empty],
                                                        r_TstarU_max[!empty],
                                                        k_min[!empty],
                                                        k_max[!empty],
                                                        obs_draws[!empty, ])
    restricted <- data.frame(beta_lower = beta_lower,
                             beta_upper = beta_upper,
                             r_uz_lower = r_uz_restricted$min,
                             r_uz_upper = r_uz_restricted$max)
    restricted <- head(restricted, -1) # Removing average row at the end
  } else {
    restricted <- NULL
  }
  unrestricted <- data.frame(k_tilde_lower = k_tilde_lower,
                             k_lower = unrestricted_bounds$k$Lower,
                             r_uz_lower = unrestricted_bounds$r_uz$Lower,
                             r_uz_upper = unrestricted_bounds$r_uz$Upper)
  list(observables = head(obs_draws, -1), #Removing extra row
       empty = head(empty, -1),
       unrestricted = head(unrestricted, -1),
       k_restriction = k_restriction,
       r_TstarU_restriction = r_TstarU_restriction,
       restricted = restricted,
       not_positive_definite = head(obs_draws$not_positive_definite, 1),
       beta_center = c(tail(beta_lower, 1), tail(beta_upper, 1)), # Adding avg beta lower and upper for interval coverage
       r_uz_center = c(tail(r_uz_restricted$min, 1), tail(r_uz_restricted$max, 1)))
}

draw_posterior <- function(y_name, T_name, z_name, data, controls = NULL,
                           r_TstarU_restriction, k_restriction = NULL,
                           n_RF_draws = 1000, n_IS_draws = 1000,
                           resample = FALSE) {

  obs_draws <- draw_observables(y_name, T_name, z_name, data, controls, n_RF_draws)
  k_tilde_lower <- get_bounds_unrest(obs_draws)$k_tilde$Lower
  k_min <- pmax(min(k_restriction), k_tilde_lower)
  k_max <- rep(min(max(k_restriction), 1), n_RF_draws)

  r_TstarU_min <- rep(min(r_TstarU_restriction), n_RF_draws)
  r_TstarU_max <- rep(max(r_TstarU_restriction), n_RF_draws)

  # The identified set is empty whenever (k_max < k_tilde_lower) in which case
  # we don't make any posterior draws
  empty <- (k_max < k_tilde_lower)
  nonempty_sets <- which(!empty)
  posterior_draws <- array(NA_real_, dim = c(n_IS_draws, 5, sum(!empty)))

  # Separate index for third dimension of posterior_draws
  posterior_draws_index <- 1
  for (i in nonempty_sets) {
    obs <- obs_draws[i, ]
    k_tilde <- runif(n_IS_draws, k_min[i], k_max[i])
    r_TstarU <- runif(n_IS_draws, r_TstarU_min[i], r_TstarU_max[i])

    if (resample) {
      M <- get_M(r_TstarU, k_tilde, obs)
      # Note: weights for sample need *not* sum to one
      random_indices <- sample(seq_len(n_IS_draws), size = n_IS_draws,
                               replace = TRUE, prob = M / max(M))
      k_tilde <- k_tilde[random_indices]
      r_TstarU <- r_TstarU[random_indices]
    }

    r_uz_final <- get_r_uz(r_TstarU, k_tilde, obs)
    s_u_final <- get_s_u(r_TstarU, k_tilde, obs)
    beta_final <- get_beta(r_TstarU, k_tilde, obs)

    posterior_draws[, , posterior_draws_index] <- cbind(
        r_TstarU,
        with(obs, (1 - T_Rsq) * k_tilde + T_Rsq), # kappa
        r_uz_final, s_u_final, beta_final)
    posterior_draws_index <- posterior_draws_index + 1
  }
  posterior_draws <- collapse_3d_array(posterior_draws)
  colnames(posterior_draws) <- c("r_TstarU", "k", "r_uz", "s_u", "beta")
  posterior_draws <- as.data.frame(posterior_draws)

  # Getting beta center
  mean_obs <- lapply(obs_draws, mean) # Adding mean obs_draws for centering
  k_tilde_lower <- get_bounds_unrest(mean_obs)$k_tilde$Lower
  k_min <- pmax(min(k_restriction), k_tilde_lower)
  k_max <- rep(min(max(k_restriction), 1), n_RF_draws)

  r_TstarU_min <- rep(min(r_TstarU_restriction), n_RF_draws)
  r_TstarU_max <- rep(max(r_TstarU_restriction), n_RF_draws)
  k_tilde <- runif(n_IS_draws, k_min, k_max)
  r_TstarU <- runif(n_IS_draws, r_TstarU_min, r_TstarU_max)

  if (resample) {
    M <- get_M(r_TstarU, k_tilde, obs)
    # Note: weights for sample need *not* sum to one
    random_indices <- sample(seq_len(n_IS_draws), size = n_IS_draws,
                             replace = TRUE, prob = M / max(M))
    k_tilde <- k_tilde[random_indices]
    r_TstarU <- r_TstarU[random_indices]
  }
  beta_centers <- get_beta(r_TstarU, k_tilde, obs)
  beta_center <- c(min(beta_centers), max(beta_centers))

  list(observables = obs_draws,
       k_restriction = k_restriction,
       r_TstarU_restriction = r_TstarU_restriction,
       posterior = posterior_draws,
       not_positive_definite = obs_draws$not_positive_definite,
       beta_center = beta_center)
}
