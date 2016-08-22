draw_sigma_CLT <- function(y, Tobs, z, n_draws) {
  n <- length(y)
  e_x <- resid(lm(y ~ Tobs))
  e_z <- resid(lm(y ~ z))
  V <- cov(cbind(Tobs * e_x, z * e_z))
  Sigma <- cov(cbind(Tobs, y, z))
  sims <- MASS::mvrnorm(n_draws, c(Sigma['Tobs', 'y'], Sigma['z', 'y']), V / n)
  g <- function(sim_row){
    out <- Sigma
    out['Tobs', 'y'] <- out['y','Tobs'] <- sim_row[1]
    out['y', 'z'] <- out['z', 'y'] <- sim_row[2]
    return(out)
  }
  Sigma_draws <- apply(sims, 1, g)
  dim(Sigma_draws) <- c(3, 3, n_draws)
  if (any(sapply(toList(Sigma_draws), det) < 0)) {
    stop("Error: non-positive definite covariance matrix drawn")
  }
  rownames(Sigma_draws) <- colnames(Sigma_draws) <- c('Tobs', 'y', 'z')
  return(Sigma_draws)
}

draw_sigma_jeffreys <- function(y, Tobs, z, n_draws) {
  n <- length(y)
  S <- (n - 1) * cov(cbind(Tobs, y, z))
  Sigma_draws <- rinvwish(n_draws, n - 1, S)
  rownames(Sigma_draws) <- colnames(Sigma_draws) <- c('Tobs', 'y', 'z')
  return(Sigma_draws)
}

draw_observables <- function(y_name, T_name, z_name, data, controls = NULL,
                             n_draws = 5000, Jeffreys = FALSE) {

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

  if(Jeffreys) {
    Sigma_draws <- draw_sigma_jeffreys(y, Tobs, z, n_draws)
  } else {
    Sigma_draws <- draw_sigma_CLT(y, Tobs, z, n_draws)
  }

  s2_T <- Sigma_draws['Tobs', 'Tobs', ]
  s2_y <- Sigma_draws['y', 'y', ]
  s2_z <- Sigma_draws['z', 'z', ]
  s_T <- sqrt(s2_T)
  s_y <- sqrt(s2_y)
  s_z <- sqrt(s2_z)
  s_Ty <- Sigma_draws['Tobs', 'y', ]
  s_Tz <- Sigma_draws['Tobs', 'z', ]
  s_zy <- Sigma_draws['z', 'y', ]
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
            r_zy = r_zy)
}

draw_bounds <- function(y_name, T_name, z_name, data, controls = NULL,
                        r_TstarU_restriction = NULL, k_restriction = NULL,
                        n_draws = 5000, Jeffreys = FALSE) {

  obs_draws <- draw_observables(y_name, T_name, z_name, data, controls,
                                n_draws, Jeffreys)

  k_tilde_lower <- get_k_tilde_lower(obs_draws)

  if (!is.null(r_TstarU_restriction)) {
    r_TstarU_min <- min(r_TstarU_restriction)
    r_TstarU_max <- max(r_TstarU_restriction)

    if (!is.null(k_restriction)) {
      # User states beliefs over kappa but we work with kappa_tilde
      T_Rsq <- obs_draws$T_Rsq[1] # All elements of obs_draws$T_Rsq are the same
      k_min <- (min(k_restriction) - T_Rsq) / (1 - T_Rsq)
      k_max <- (max(k_restriction) - T_Rsq) / (1 - T_Rsq)
      k_min <- pmax(k_tilde_lower, k_min) # vector: could vary with obs_draws row
      k_max <- min(1, k_max) # always a scalar
    } else {
      k_min <- k_tilde_lower # vector: varies with obs_draws row
      k_max <- 1 # always a scalar
    }

    # We ensure above that but the user may have specified a k_max that is less
    # than some elements of k_tilde_lower as in one of the examples for Colonial
    # Origins from the paper. When this occurs, the identified set is empty.
    empty <- k_max < k_tilde_lower

    # Only compute the bounds for the non-empty identified sets
    beta_lower <- get_beta_lower(r_TstarU_max, k_min[!empty], k_max,
                                 obs_draws[!empty, ])
    beta_upper <- get_beta_upper(r_TstarU_min, k_min[!empty], k_max,
                                 obs_draws[!empty, ])
    #rho_uz_lower_restricted <- # Not yet implemented
    #rho_uz_upper_restricted <- # Not yet implemented
    restricted <- data.frame(beta_lower = beta_lower,
                             beta_upper = beta_upper)
                             #r_uz_lower = r_uz_lower_restricted,
                             #r_uz_upper = r_uz_upper_restricted)
  } else {
    restricted <- NULL
  }
  unrestricted = data.frame(k_tilde_lower = k_tilde_lower,
                            k_lower = get_k_lower(obs_draws),
                            r_uz_lower = get_r_uz_lower(obs_draws),
                            r_uz_upper = get_r_uz_upper(obs_draws))
  list(observables = obs_draws,
       empty = empty,
       unrestricted = unrestricted,
       k_restriction = k_restriction,
       r_TstarU_restriction = r_TstarU_restriction,
       restricted = restricted)
}


draw_posterior <- function(y_name, T_name, z_name, data, controls = NULL,
                           r_TstarU_restriction, k_restriction = NULL,
                           n_RF_draws = 1000, n_IS_draws = 1000,
                           Jeffreys = FALSE, resample = FALSE) {

  obs_draws <- draw_observables(y_name, T_name, z_name, data, controls,
                                n_RF_draws, Jeffreys)

  k_tilde_lower <- get_k_tilde_lower(obs_draws)
  k_min <- pmax(min(k_restriction), k_tilde_lower)
  k_max <- min(max(k_restriction), 1)

  r_TstarU_min <- min(r_TstarU_restriction)
  r_TstarU_max <- max(r_TstarU_restriction)

  # The identified set is empty whenever (k_max < k_tilde_lower) in which case
  # we don't make any posterior draws
  empty <- k_max < k_tilde_lower
  nonempty_sets <- which(!empty)
  posterior_draws <- array(NA_real_, dim = c(n_IS_draws, 5, sum(!empty)))

  # Separate index for third dimension of posterior_draws
  posterior_draws_index <- 1
  for (i in nonempty_sets) {

    obs <- obs_draws[i, ]
    k_tilde <- runif(n_IS_draws, k_min[i], k_max)
    r_TstarU <- runif(n_IS_draws, r_TstarU_min, r_TstarU_max)

    if (resample) {
      M <- get_M(r_TstarU, k_tilde, obs)
      # Note: weights for sample need *not* sum to one
      random_indices <- sample(seq_len(n_IS_draws), size = n_IS_draws,
                               replace = TRUE, prob = M / max(M))
      k_tilde <- k_tilde[random_indices]
      r_TstarU <- r_TstarU[random_indices]
    }

    posterior_draws[, , posterior_draws_index] <- cbind(
        r_TstarU,
        with(obs, (1 - T_Rsq[1]) * k_tilde + T_Rsq[1]), # kappa
        get_r_uz(r_TstarU, k_tilde, obs),
        get_s_u(r_TstarU, k_tilde, obs),
        get_beta(r_TstarU, k_tilde, obs))
    posterior_draws_index <- posterior_draws_index + 1
  }
  posterior_draws <- collapse_3d_array(posterior_draws)
  colnames(posterior_draws) <- c('r_TstarU', 'k', 'r_uz', 's_u', 'beta')
  posterior_draws <- as.data.frame(posterior_draws)

  list(observables = obs_draws,
       empty = empty,
       k_restriction = k_restriction,
       r_TstarU_restriction = r_TstarU_restriction,
       posterior = posterior_draws)
}
