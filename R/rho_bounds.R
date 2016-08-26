r_TstarU_upper <- 0.9
r_TstarU_lower <- 0
k_lower <- rep(0.7, nrow(obs))
k_upper <- k_upper <- rep(1, nrow(obs))

# Candidate Set I - Corner for both kappa and rho_TstarU
corner1 <- get_r_uz(r_TstarU_upper, k_upper, obs)
corner2 <- get_r_uz(r_TstarU_lower, k_lower, obs)
corner3 <- get_r_uz(r_TstarU_upper, k_lower, obs)
corner4 <- get_r_uz(r_TstarU_lower, k_upper, obs)
min_corner <- pmin(corner1, corner2, corner3, corner4)
max_corner <- pmax(corner1, corner2, corner3, corner4)

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

min_corner_r <- pmin(min_corner_r1, min_corner_r2)
max_corner_r <- pmax(max_corner_r1, max_corner_r2)

# Finally: overall max and min
r_uz_max <- pmax(max_corner, max_corner_r, max_corner_k)
r_uz_min <- pmin(min_corner, min_corner_r, min_corner_k)
