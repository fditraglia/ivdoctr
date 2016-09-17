get_observables <- function(y_name, T_name, z_name, data, controls = NULL) {

  # Project out control regressors if present
  if (!is.null(controls)) {
    y <- resid(lm(reformulate(controls, response = y_name), data))
    T_reg <- lm(reformulate(controls, response = T_name), data)
    z <- resid(lm(reformulate(controls, response = z_name), data))
    Tobs <- resid(T_reg)
    T_Rsq <- summary(T_reg)$r.squared
    rm(T_reg, z_reg)
  } else {
    y <- get(y_name, data)
    Tobs <- get(T_name, data)
    z <- get(z_name, data)
    T_Rsq <- 0 # No controls is the same as controls that are uncorrelated with Tobs
  }

  Sigma <- cov(cbind(Tobs, y, z))
  s2_T <- Sigma['Tobs', 'Tobs']
  s2_y <- Sigma['y', 'y']
  s2_z <- Sigma['z', 'z']
  s_Ty <- Sigma['Tobs', 'y']
  s_Tz <- Sigma['Tobs', 'z']
  s_zy <- Sigma['z', 'y']

  Rho <- cov2cor(Sigma)
  r_Ty <- Rho['Tobs', 'y']
  r_Tz <- Rho['Tobs', 'z']
  r_zy <- Rho['z', 'y']

  ans <- list(n = nrow(data),
              T_Rsq = T_Rsq,
              s2_T = s2_T,
              s2_y = s2_y,
              s2_z = s2_z,
              s_Ty = s_Ty,
              s_Tz = s_Tz,
              s_zy = s_zy,
              r_Ty = r_Ty,
              r_Tz = r_Tz,
              r_zy = r_zy)
  return(ans)
}

get_k_tilde_lower <- function(obs) {
  ans <- with(obs, (r_Ty^2 + r_Tz^2 - 2 * r_Ty * r_Tz * r_zy) / (1 - r_zy^2))
  return(ans)
}

get_k_lower <- function(obs) {
  k_tilde_lower <- get_k_tilde_lower(obs)
  ans <- with(obs, (1 - T_Rsq) * k_tilde_lower + T_Rsq)
  return(ans)
}

get_r_uz_lower <- function(obs) {
  k_tilde_lower <- get_k_tilde_lower(obs)
  nontrivial_lower_bound <- with(obs, r_Ty * r_Tz < k_tilde_lower * r_zy)
  ifelse(nontrivial_lower_bound, -1 * abs(obs$r_Tz) / sqrt(k_tilde_lower), -1)
}

get_r_uz_upper <- function(obs) {
  k_tilde_lower <- get_k_tilde_lower(obs)
  nontrivial_upper_bound <- with(obs, r_Ty * r_Tz >= k_tilde_lower * r_zy)
  ifelse(nontrivial_upper_bound, abs(obs$r_Tz) / sqrt(k_tilde_lower), 1)
}

get_r_uz <- function(r_TstarU, k, obs) {
  A <- with(obs, r_TstarU * r_Tz / sqrt(k))
  B1 <- with(obs, r_Ty * r_Tz - k * r_zy)
  B2 <- with(obs, sqrt((1 - r_TstarU^2) / (k * (k - r_Ty^2))))
  A - B1 * B2
}

get_M <- function(r_TstarU, k, obs) {
  A <- with(obs, r_TstarU * r_Tz / sqrt(k))
  B1 <- with(obs, r_Ty * r_Tz - k * r_zy)
  B2 <- with(obs, sqrt((1 - r_TstarU^2) / (k * (k - r_Ty^2))))
  dr_TstarU <- with(obs, r_Tz / sqrt(k) + r_TstarU * B1 /
                    sqrt(k * (k - r_Ty^2) * (1 - r_TstarU^2)))
  dk <- with(obs, -r_TstarU * r_Tz / (2 * k^(3 / 2)) + B2 *
                  (r_zy + B1 / 2.0 * (1 / k + 1/(k - r_Ty^2))))
  sqrt(1 + dr_TstarU^2 + dk^2)
}

get_s_u <- function(r_TstarU, k, obs) {
  with(obs, sqrt(s2_y * (k - r_Ty^2) / (k * (1 - r_TstarU^2))))
}

get_beta <- function(r_TstarU, k, obs) {
  r_uz <- get_r_uz(r_TstarU, k, obs)
  s_u <- get_s_u(r_TstarU, k, obs)
  with(obs, (r_zy * sqrt(s2_y) - r_uz * s_u) / (r_Tz * sqrt(s2_T)))
}

# This function is vectorized wrt k_min and and obs (k_max is always a scalar).
get_beta_lower <- function(r_TstarU_max, k_min, k_max, obs) {
  # min(beta) occurs at max(r_TstarU) but could be a corner value for kappa
  beta_corner <- pmin(get_beta(r_TstarU_max, k_min, obs),
                      get_beta(r_TstarU_max, k_max, obs))

  # If r_TstarU_max = 0, no need to check the interior solution
  if (identical(r_TstarU_max, 0)) {
    out <- beta_corner
  } else {
    C <- 1 / sqrt(1 + (r_TstarU_max^2 / (1 - r_TstarU_max^2)))
    k_interior <- with(obs, 2 * r_Ty^2 / (1 - sign(r_TstarU_max * r_Ty) * C))

    # It only makes sense to calculate beta_interior if k_interior is between
    # k_min and k_max. If not, set it to Inf so it always exceeds beta_corner.
    k_interior_is_valid <- (k_interior > k_min) & (k_interior < k_max)
    beta_interior <- ifelse(k_interior_is_valid,
                            get_beta(r_TstarU_max, k_interior, obs), Inf)
    out <- pmin(beta_corner, beta_interior)
  }
  return(out)
}

# This function is vectorized wrt k_min and and obs (k_max is always a scalar).
get_beta_upper <- function(r_TstarU_min, k_min, k_max, obs) {
  # max(beta) occurs at min(r_TstarU) but could be a corner value for kappa
  beta_corner <- pmax(get_beta(r_TstarU_min, k_min, obs),
                      get_beta(r_TstarU_min, k_max, obs))

  # If r_TstarU_min = 0, no need to check interior solution
  if(identical(r_TstarU_min, 0)) {
    out <- beta_corner
  } else {
    C <- 1 / sqrt(1 + (r_TstarU_min^2 / (1 - r_TstarU_min^2)))
    k_interior <- with(obs, 2 * r_Ty^2 / (1 + sign(r_TstarU_min * r_Ty) * C))

    # It only makes sense to calculate beta_interior if k_interior is between
    # k_min and k_max. If not, set it to -Inf so it never exceeds beta_corner.
    k_interior_is_valid <- (k_interior > k_min) & (k_interior < k_max)
    beta_interior <- ifelse(k_interior_is_valid,
                            get_beta(r_TstarU_min, k_interior, obs), -Inf)
    out <- pmax(beta_corner, beta_interior)
  }
  return(out)
}
