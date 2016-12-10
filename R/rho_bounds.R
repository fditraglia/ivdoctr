get_r_uz_bounds <- function(r_TstarU_lower, r_TstarU_upper,
                            k_lower, k_upper, obs) {
  # If beliefs are not restrictive, function returns unrestricted bounds
  unrestricted <- get_bounds_unrest(obs)
  if (r_TstarU_lower <= unrestricted$r_TstarU$Lower &
      r_TstarU_upper >= unrestricted$r_TstarU$Upper &
      k_lower <= unrestricted$k$Lower &
      k_upper >= unrestricted$k$Upper) {
    warning("Beliefs are not restrictive, returning unrestricted bounds")
    return(unrestricted)
  } else {

    # Storing the more restrictive bounds for each parameter 
    # (either user or unrestricted)
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
    
    # Finding r_uz bounds based on restricted bounds above
    # Candidate Set I - Corner for both kappa and rho_TstarU
    set1 <- candidate1(r_TstarU_lower_bound, r_TstarU_upper_bound, 
                       k_lower_bound, k_upper_bound, obs)
  
    # Candidate Set II - Corner for kappa, interior for rho_TstarU
    set2 <- candidate2(r_TstarU_lower_bound, r_TstarU_upper_bound, 
                       k_lower_bound, k_upper_bound, obs)

    # Candidate Set III - Interior for kappa, corner for rho_TstarU
    set3 <- candidate3(r_TstarU_lower_bound, r_TstarU_upper_bound,
                       k_lower_bound, k_upper_bound, obs)
    
    # Finally: overall max and min
    r_uz_max <- pmax(set1$max_corner, set2$max_edge, set3$r_uz$max_edge, na.rm = TRUE)
    r_uz_min <- pmin(set1$min_corner, set2$min_edge, set3$r_uz$min_edge, na.rm = TRUE)
    data.frame(min = r_uz_min, max = r_uz_max)
  }
}

# Evaluates the corners given user bounds. Vectorized wrt r_TstarU and k bounds.
candidate1 <- function(r_TstarU_lower, r_TstarU_upper, k_lower, k_upper, obs) {
    corner1 <- get_r_uz(r_TstarU_upper, k_upper, obs)
    corner2 <- get_r_uz(r_TstarU_lower, k_lower, obs)
    corner3 <- get_r_uz(r_TstarU_upper, k_lower, obs)
    corner4 <- get_r_uz(r_TstarU_lower, k_upper, obs)
    min_corner <- pmin(corner1, corner2, corner3, corner4, na.rm = TRUE)
    max_corner <- pmax(corner1, corner2, corner3, corner4, na.rm = TRUE)
    ans <- list(min_corner = min_corner, 
                max_corner = max_corner)
    return(ans)
}

# Evaluates the edge where k is on the boundary. Vectorized wrt r_TstarU and k bounds.
candidate2 <- function(r_TstarU_lower, r_TstarU_upper, k_lower, k_upper, obs) {
    k_bounds <- cbind(k_lower, k_upper)
    a <- with(obs, r_Tz * sqrt(k_bounds - r_Ty^2) / (r_Ty * r_Tz - k_bounds * r_zy))
    r <- -1 * a / sqrt(1 + a^2)
    r <- ifelse(r <= r_TstarU_upper & r >= r_TstarU_lower, r, NA)
    edges <- get_r_uz(r, k_bounds, obs)
    min_edge <- pmin(edges[, 1], edges[, 2], na.rm = TRUE)
    max_edge <- pmax(edges[, 1], edges[, 2], na.rm = TRUE)
    ans <- list(min_edge = min_edge, 
                max_edge = max_edge)
    return(ans)
}

# Evaluates the edge where r_TstarU is on the boundary. 
candidate3 <- function(r_TstarU_lower, r_TstarU_upper, k_lower, k_upper, obs) {
  r_TstarU_bounds <- c(r_TstarU_lower, r_TstarU_upper)
  # Getting cubic coefficients
  d <- with(obs, -r_Tz ^ 2 * r_TstarU_bounds ^ 2)
  c <- with(obs, 3 * r_Ty ^ 2 * r_Tz ^ 2 * r_TstarU_bounds ^ 2 + 
                     ((2 * r_Ty * r_Tz - r_zy * r_Ty ^ 2) ^ 2) * 
                     (1 - r_TstarU_bounds ^ 2))
  b <- with(obs, -3 * r_Ty ^ 4 * r_Tz ^ 2 * r_TstarU_bounds ^ 2 - 
                  2 * (2 * r_Ty * r_Tz - r_zy * r_Ty ^ 2) * 
                      (r_Tz * r_Ty ^ 3) * (1 - r_TstarU_bounds ^ 2))
  a <- rep(with(obs, r_Ty ^ 6 * r_Tz ^ 2), 2)
  coefs <- cbind(a, b, c, d)
  
  # Getting roots of cubic
  all_roots <- apply(coefs, 1, polyroot)
  real_roots <- Re(unlist(all_roots))
  real_roots <- real_roots[(real_roots >= k_lower) & (real_roots <= k_upper)]

  # Evaluating all combinations of roots and bounds for r_TstarU
  r_uz <- outer(real_roots, r_TstarU_bounds, 
                function(real_roots, r_TstarU_bounds) 
                  get_r_uz(r_TstarU_bounds, real_roots, obs))
  
  # Getting max and min if there are real roots
  min_edge <- ifelse(length(real_roots) > 0, min(r_uz), NA)
  max_edge <- ifelse(length(real_roots) > 0, max(r_uz), NA)
  ans <- list(r_uz = list(min_edge = min_edge, max_edge = max_edge),
              k_roots = real_roots)
  return(ans)
}
