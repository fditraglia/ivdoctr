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
  Sigma_draws <- toList(Sigma_draws)
  if (any(sapply(Sigma_draws, det) < 0)) {
    stop("Error: non-positive definite covariance matrix drawn")
  }
  return(Sigma_draws)
}

draw_sigma_jeffreys <- function(y, Tobs, z, n_draws) {
  n <- length(y)
  S <- (n - 1) * cov(cbind(Tobs, y, z))
  Sigma_draws <- rinvwish(n_draws, n - 1, S)
  toList(Sigma_draws)
}

draw_observables <- function(y_name, T_name, z_name, controls = NULL, data,
                             nDraws = 5000, Jeffreys = FALSE) {

}
