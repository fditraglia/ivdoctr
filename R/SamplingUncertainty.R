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

draw_bounds <- function(r_TstarU_range, k_range = NULL,
                        y_name, T_name, z_name, data, controls = NULL,
                        n_draws = 5000, Jeffreys = FALSE) {

  obs_draws <- draw_observables(y_name, T_name, z_name, data, controls,
                                n_draws, Jeffreys)

  r_TstarU_min <- min(r_TstarU_range)
  r_TstarU_max <- max(r_TstarU_range)

  k_tilde_lower <- get_k_tilde_lower(obs_draws)
  k_lower <- get_k_lower(obs_draws)
  r_uz_lower <- get_r_uz_lower(obs_draws)
  r_uz_upper <- get_r_uz_upper(obs_draws)

  if (!is.null(k_range)) {
    k_min <- min(k_range)
    k_max <- max(k_range)
  } else {

  }

  # Add calculation of tighter bounds for rho_uz that incorporate beliefs later
}



