get_r_uz_bounds <- function(r_TstarU_lower, r_TstarU_upper,
                            k_lower, k_upper, obs) {
  
  # If beliefs are not restrictive, function returns unrestricted bounds
  unrestricted <- get_bounds_unrest(obs)
  if (rTstarU_lower <= unrestricted$r_TstarU$Lower &
      rTstarU_upper >= unrestricted$r_TstarU$Upper &
      k_lower <= unrestricted$k$Lower &
      k_upper >= unrestricted$k$Upper) {
    warning("Beliefs are not restrictive, returning unrestricted bounds")
    return(unrestricted)
  } else {
    
    # Storing restricted bounds for later
    r_TstarU_lower_bound <- ifelse(r_TstarU_lower <= unrestricted$r_TstarU$Lower,
                                   unrestricted$r_TstarU$Lower,
                                   r_TstarU_lower)
    r_TstarU_upper_bound <- ifelse(r_TstarU_upper >= unrestricted$r_TstarU$Upper,
                                   unrestricted$r_TstarU$Upper,
                                   r_TstarU_upper)
    k_lower_bound <- ifelse(k_lower <= unrestricted$k$Lower, 
                            unrestricted$k$Lower, 
                            k_lower)
    k_upper_bound <- ifelse(k_upper >= unrestricted$k$Upper, 
                            unrestricted$k$Upper, 
                            k_upper)
    
    # Finding r_uz bounds
    # Candidate Set I - Corner for both kappa and rho_TstarU
    set1 <- candidate1(r_TstarU_upper, r_TstarU_lower, k_upper, k_lower, obs)
  
    # Candidate Set II - Corner for kappa, interior for rho_TstarU
    a1 <- with(obs, r_Tz * sqrt(k_lower - r_Ty^2) / (r_Ty * r_Tz - k_lower * r_zy))
    r1 <- -1 * a1 / sqrt(1 + a1^2)
    r1_ok <- (r1 > r_TstarU_lower) & (r1 < r_TstarU_upper)
    corner_k1 <- rep(NA_real_, length(r1))
    corner_k1[r1_ok] <- get_r_uz(r1[r1_ok], k_lower[r1_ok], obs[r1_ok,])
  
    a2 <- with(obs, r_Tz * sqrt(k_upper - r_Ty^2) / (r_Ty * r_Tz - k_upper * r_zy))
    r2 <- -1 * a2 / sqrt(1 + a2^2)
    r2_ok <- (r2 > r_TstarU_lower) & (r2 < r_TstarU_upper)
    corner_k2 <- rep(NA_real_, length(r2))
    corner_k2[r2_ok] <- get_r_uz(r2[r2_ok], k_lower[r2_ok], obs[r2_ok,])
  
    min_corner_k <- pmin(corner_k1, corner_k2, na.rm = TRUE)
    max_corner_k <- pmax(corner_k1, corner_k2, na.rm = TRUE)
  
    # Candidate Set III - Interior for kappa, corner for rho_TstarU
    c0 <- with(obs, -0.25 * r_Ty^6 * r_Tz^2)
    c1_1 <- with(obs, r_Ty^5 * r_Tz * r_zy)
    c1_2 <- with(obs, r_Ty^4 * r_Tz^2)
    c2_1 <- with(obs, r_Ty^4 * r_zy^2)
    c2_2 <- with(obs, r_Ty^3 * r_Tz * r_zy)
    c2_3 <- with(obs, r_Ty^2 * r_Tz^2)
  
    c1_lower <- -0.5 * c1_1 + c1_2 + 0.25 * (2 * c1_1 - c1_2) * r_TstarU_lower
    c2_lower <- -0.25 * c2_1 + c2_2 - c2_3 +
      0.25 * (c2_1 - 4 * c2_2 + c2_3) * r_TstarU_lower
    c3_lower <- with(obs, 0.25 * r_TstarU_lower^2 * r_Tz^2)
    coefs_lower <- cbind(c0, c1_lower, c2_lower, c3_lower)
    k1 <- Re(t(apply(coefs_lower, 1, polyroot)))
    k1_ok <- apply(k1, 2, function(x) (x < k_upper) & (x > k_lower))
    k1[!k1_ok] <- NA
    corner_r1 <- apply(k1, 2, function(x) get_r_uz(r_TstarU_lower, x, obs))
    min_corner_r1 <- suppressWarnings(apply(corner_r1, 1, min, na.rm = TRUE))
    max_corner_r1 <- suppressWarnings(apply(corner_r1, 1, max, na.rm = TRUE))
  
    c1_upper <- -0.5 * c1_1 + c1_2 + 0.25 * (2 * c1_1 - c1_2) * r_TstarU_upper
    c2_upper <- -0.25 * c2_1 + c2_2 - c2_3 +
      0.25 * (c2_1 - 4 * c2_2 + c2_3) * r_TstarU_upper
    c3_upper <- with(obs, 0.25 * r_TstarU_upper^2 * r_Tz^2)
    coefs_upper <- cbind(c0, c1_upper, c2_upper, c3_upper)
    k2 <- Re(t(apply(coefs_upper, 1, polyroot)))
    k2_ok <- apply(k2, 2, function(x) (x < k_upper) & (x > k_lower))
    k2[!k2_ok] <- NA
    corner_r2 <- apply(k2, 2, function(x) get_r_uz(r_TstarU_upper, x, obs))
    min_corner_r2 <- suppressWarnings(apply(corner_r2, 1, min, na.rm = TRUE))
    max_corner_r2 <- suppressWarnings(apply(corner_r2, 1, max, na.rm = TRUE))
  
    min_corner_r <- pmin(min_corner_r1, min_corner_r2, na.rm = TRUE)
    max_corner_r <- pmax(max_corner_r1, max_corner_r2, na.rm = TRUE)
  
    # Finally: overall max and min
    r_uz_max <- pmax(max_corner, max_corner_r, max_corner_k, na.rm = TRUE)
    r_uz_min <- pmin(min_corner, min_corner_r, min_corner_k, na.rm = TRUE)
    data.frame(min = r_uz_min, max = r_uz_max)
  }
}

# Evaluates the corners given user bounds. Vectorized wrt r_TstarU and k bounds.
candidate1 <- function(r_TstarU_upper, r_TstarU_lower, k_upper, k_lower, obs) {
    corner1 <- get_r_uz(r_TstarU_upper, k_upper, obs)
    corner2 <- get_r_uz(r_TstarU_lower, k_lower, obs)
    corner3 <- get_r_uz(r_TstarU_upper, k_lower, obs)
    corner4 <- get_r_uz(r_TstarU_lower, k_upper, obs)
    min_corner <- pmin(corner1, corner2, corner3, corner4, na.rm = TRUE)
    max_corner <- pmax(corner1, corner2, corner3, corner4, na.rm = TRUE)
    ans <- cbind(min_corner, max_corner)
    return(ans)
}

# Evaluates the edge where k is on the boundary. Vectorized wrt k bounds.
candidate2 <- function(r_TstarU_upper, r_TstarU_lower, k_upper, k_lower, obs) {
    k_bounds <- cbind(k_lower, k_upper)
    a <- with(obs, r_Tz * sqrt(k_bounds - r_Ty^2) / (r_Ty * r_Tz - k_bounds * r_zy))
    r <- -1 * a / sqrt(1 + a^2)
    r <- ifelse(r <= r_TstarU_upper & r >= r_TstarU_lower, r, NA)
    edges <- get_r_uz(r, k_bounds, obs)
    min_edge <- pmin(edges[, 1], edges[, 2], na.rm = TRUE)
    max_edge <- pmax(edges[, 1], edges[, 2], na.rm = TRUE)
    ans <- cbind(min_edge, max_edge)
    return(ans)
}
