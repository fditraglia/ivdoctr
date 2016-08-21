get_observables <- function(y_name, T_name, z_name, controls = NULL, data) {

  first_stage <- reformulate(c(z_name, controls), response = NULL)
  second_stage <- reformulate(c(T_name, controls), response = y_name)
  OLS <- lm(second_stage, data)
  IV <-  AER::ivreg(second_stage, first_stage, data)

  b_OLS <- coef(OLS)[T_name]
  b_IV <- coef(IV)[T_name]
  se_OLS <- sqrt(diag(vcov(OLS)))[T_name]
  se_IV <- sqrt(diag(vcov(IV)))[T_name]

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

  Sigma <- cov(cbind(Tobs, y, z))
  s2_T <- Sigma['Tobs', 'Tobs']
  s2_y <- Sigma['y', 'y']
  s2_z <- Sigma['z', 'z']
  s_Ty <- Sigma['Tobs', 'y']
  s_Tz <- Sigma['Tobs', 'z']
  s_zy <- Sigma["z", "y"]

  Rho <- cov2cor(Sigma)
  r_Ty <- Rho['Tobs', 'y']
  r_Tz <- Rho['Tobs', 'z']
  r_zy <- Rho['z', 'y']

  k_tilde_lower <- (r_Ty^2 + r_Tz^2 - 2 * r_Ty * r_Tz * r_zy) / (1 - r_zy^2)
  k_lower <- (1 - T_Rsq) * k_tilde_lower + T_Rsq

  if (r_Ty * r_Tz < k_tilde_lower * r_zy) {
    r_uz_lower <- -1 * abs(r_Tz) / sqrt(k_tilde_lower)
    r_uz_upper <- 1
  } else {
    r_uz_lower <- -1
    r_uz_upper <- abs(r_Tz) / sqrt(k_tilde_lower)
  }

  list(n = nrow(data),
       b_OLS = b_OLS,
       se_OLS = se_OLS,
       b_IV = b_IV,
       se_IV = se_IV,
       T_Rsq = T_Rsq,
       z_Rsq = z_Rsq,
       Sigma = Sigma,
       s2_T = s2_T,
       s2_y = s2_y,
       s2_z = s2_z,
       s_Ty = s_Ty,
       s_Tz = s_Tz,
       s_zy = s_zy,
       Rho = Rho,
       r_Ty = r_Ty,
       r_Tz = r_Tz,
       r_zy = r_zy,
       k_lower = k_lower,
       k_tilde_lower = k_tilde_lower,
       r_uz_lower = r_uz_lower,
       r_uz_upper = r_uz_upper)
}

get_r_zu <- function(r_TstarU, k, obs) {
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
  r_zu <- get_r_zu(r_TstarU, k, obs)
  s_u <- get_s_u(r_TstarU, k, obs)
  with(obs, (r_zy * sqrt(s2_y) - r_zu * s_u) / (r_Tz * sqrt(s2_T)))
}


get_beta_lower <- function(r_TstarU_range, k_range, obs) {
  # Min for beta always occurs at *max* for r_TstarU
  r_TstarU_max <- max(r_TstarU_range)
  k_min <- min(k_range)
  k_max <- max(k_range)

  # Solution could be at a corner value for kappa
  beta_corner <- min(get_beta(r_TstarU_max, k_min, obs),
                     get_beta(r_TstarU_max, k_max, obs))

  # If r_TstarU_max = 0, no need to check the interior solution
  if (identical(r_TstarU_max, 0)) {
    out <- beta_corner
  } else {
    C <- 1 / sqrt(1 + (r_TstarU_max^2 / (1 - r_TstarU_max^2)))
    k_interior <- with(obs, 2 * r_Ty^2 / (1 - sign(r_TstarU_max * r_Ty) * C))

    # Only need to check beta_interior if k_interior is between k_min and k_max
    if ((k_interior > k_min) && (k_interior < k_max)) {
      beta_interior <- get_beta(r_TstarU_max, k_interior, obs)
      out <- min(beta_corner, beta_interior)
    } else {
      out <- beta_corner
    }
  }
  return(out)
}

get_beta_upper <- function(r_TstarU_range, k_range, obs) {
  # Max for beta always occurs at *min* for r_TstarU but could be interior or
  r_TstarU_min <- min(r_TstarU_range)
  k_min <- min(k_range)
  k_max <- max(k_range)

  # Solution could be at a corner for kappa
  beta_corner <- max(get_beta(r_TstarU_min, k_min, obs),
                     get_beta(r_TstarU_min, k_max, obs))

  # If r_TstarU_min = 0, no need to check interior solution
  if(identical(r_TstarU_min, 0)) {
    out <- beta_corner
  } else {
    C <- 1 / sqrt(1 + (r_TstarU_min^2 / (1 - r_TstarU_min^2)))
    k_interior <- with(obs, 2 * r_Ty^2 / (1 + sign(r_TstarU_min * r_Ty) * C))

    # Only need to check beta_interior if k_interior is between k_min and k_max
    if ((k_interior > k_min) && (k_interior < k_max)){
      beta_interior <- get_beta(r_TstarU_min, k_interior, obs)
      out <- max(beta_corner, beta_interior)
    } else {
      out <- beta_corner
    }
  }
  return(out)
}












