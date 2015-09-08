#' Convert 3-d array to list of matrixes
#'
#' @param myArray A three-dimensional numeric array.
#' @return A list of numeric matrices.
#' @examples
#' M <- array(c(1, 1, 1, 1, 2, 2, 2, 2), c(2, 2, 2))
#' toList(M)
toList <- function(myArray){
  lapply(seq_len(dim(myArray)[3]), function(i) myArray[,,i,drop = TRUE])
}

SuTilde2Su <- function(S, SuTildeVec){
#-------------------------------------------------------------
# Convert Sigma_u_tilde to Sigma_u
#-------------------------------------------------------------
  Sy <- sqrt(S[2,2])
  return(SuTildeVec * Sy)
}

getBetaIV <- function(S){
#-------------------------------------------------------------
# Given cov matrix S[x,y,z] calculate IV estimator
#-------------------------------------------------------------
  Sxz <- S[1,3]
  Syz <- S[2,3]
  return(Syz / Sxz)
}

getBetaOLS <- function(S){
#-------------------------------------------------------------
# Given cov matrix S[x,y,z] calculate OLS estimator
#-------------------------------------------------------------
  Sxy <- S[1,2]
  Sx2 <- S[1,1]
  return(Sxy / Sx2)
}

getBeta <- function(S, SuVec, RzuVec){
#-------------------------------------------------------------
# Given cov matrix S[x,y,z] and values for Sigma_u, Rho_zu
# calculate the implied "true" Beta
#-------------------------------------------------------------
  BetaIV <- getBetaIV(S)
  Sxz <- S[1,3]
  Sz2 <- S[3,3]
  PI <- Sxz / Sz2
  numerator <- SuVec * RzuVec
  denominator <- PI * sqrt(Sz2)
  return(BetaIV - numerator / denominator)
}

# Convert from Ktilde to K
toKappa <- function(Ktilde, xRsq){
  Ktilde * (1 - xRsq) + xRsq
}

# Convert from RzuTilde to Rzu
toRzu <- function(RzuTilde, zRsq){
   RzuTilde * sqrt(1 - zRsq)
 }

# Convert from RxsuTilde to Rxsu
toRxsu <- function(RxsuTilde, xRsq, Kappa){
   RxsuTilde * sqrt(1 - xRsq / Kappa)
}

#' Posterior draws for a cov matrix based on Jeffrey's prior
#'
#' @param inData A matrix or dataframe, each column of which is a variable and
#' each row an observation.
#' @param n_draws An integer, the number of posterior draws.
#' @return An array of numeric matrices, each of which is a sample from the
#' marginal posterior of the covariance matrix.
#' @details Based on a joint normal likelihood, not conditional as in
#' regression.
#' @examples
#' set.seed(1234)
#' x <- rnorm(100)
#' y <- x / 3 + rnorm(100)
#' z <- y / 5 + rnorm(100)
#' covJeffreys(cbind(x, y, z), 10)
covJeffreys <- function(inData, n_draws){
  n <- nrow(inData)
  S <- (n - 1) * cov(inData)
  return(rinvwish(n_draws, n - 1, S))
}


#' Posterior draws for IV and OLS based on a large-sample approximation
#'
#' @param inData A dataframe with columns named x, y, and z each of which is a
#' numeric vector.
#' @param n_draws An integer, the number of posterior draws.
#' @return An array of numeric matrices, each of which is a sample from the
#' marginal posterior of the covariance matrix.
#' @details In OLS regression we condition on x and in IV regression we
#' condition on x and z. This function generates posterior draws for the
#' covariance matrix of (x, y, z) conditional on the observed values of x and
#' z based on the usual large-sample frequentist asymptotics. In particular
#' from the CLT, the OLS and IV estimators are jointly asymptotically normal.
#' Treating this asymptotic result as exact and multiplying through by the
#' sample estimates of Var(x) and Cov(x,z) gives a joint normal distribution for
#' Cov(x,y) and Cov(z,y). To get draws for Cov(x,y,z) we simply fill in the
#' sample estimates of Var(x), Var(y), Var(z), and Cov(x,z). When generating
#' draws for Cov(x,y) and Cov(z,y) we allow for heteroskedasticity. Note that
#' when generated in this way, posterior draws are not guaranteed to be positive
#' definite. If any of the draws is not positive definite the function throws
#' an error.
#' @examples
#' set.seed(1234)
#' x <- rnorm(100)
#' y <- x / 3 + rnorm(100)
#' z <- y / 5 + rnorm(100)
#' covCLT(data.frame(x, y, z), 10)
covCLT <- function(inData, n_draws){
  stopifnot(setequal(names(inData), c("x", "y", "z")))
  x <- inData$x
  y <- inData$y
  z <- inData$z
  e_x <- lm(y ~ x)$residuals
  e_z <- lm(y ~ z)$residuals
  V <- cov(cbind(x * e_x, z * e_z))
  Sigma <- cov(inData)
  sims <- MASS::mvrnorm(n_draws, c(Sigma[1,2], Sigma[3,2]), V / length(e_z))
  upper_det <- Sigma[1,1] * Sigma[2,2] - sims[,1]^2
  if(any(upper_det <= 0)) stop("non-positive definite cov matrix generated (1)")
  full_det <- cov(x,z) * (sims[,1] * sims[,2] - cov(x,z) * var(y)) -
              sims[,2] * (var(x) * sims[,2] - cov(x,z) * sims[,1]) +
              var(z) * upper_det
  if(any(full_det <= 0)) stop("non-positive definite cov matrix generated (2)")
  g <- function(sim_row){
    out <- Sigma
    out[1,2] <- out[2,1] <- sim_row[1]
    out[2,3] <- out[3,2] <- sim_row[2]
    return(out)
  }
  Sigma_draws <- apply(sims, 1, g)
  dim(Sigma_draws) <- c(3, 3, n_draws)
  return(Sigma_draws)
}

#======================= Helper functions for plotting

get_Rzu <- function(Sigma, Rxsu, K){
  R <- cov2cor(Sigma)
  Rxy <- R[1,2]
  Rxz <- R[1,3]
  Rzy <- R[2,3]
  A <- Rxsu * Rxz / sqrt(K)
  B1 <- (Rxy * Rxz - K * Rzy)
  B2 <- sqrt((1 - Rxsu^2) / (K * (K - Rxy^2)))
  return(A - B1 * B2)
}

facet_val <- function(m){
  (m[-1,-1] + m[-1, -ncol(m)] + m[-nrow(m),-1] + m[-nrow(m), -ncol(m)])/4
}

get_Su <- function(Sigma, Rxsu, K){
  R <- cov2cor(Sigma)
  Rxy <- R[1,2]
  Rxz <- R[1,3]
  Rzy <- R[2,3]
  Sy <- sqrt(Sigma[2,2])
  A <- Rxsu * Rxz / sqrt(K)
  B1 <- (Rxy * Rxz - K * Rzy)
  B2 <- sqrt((1 - Rxsu^2) / (K * (K - Rxy^2)))
  Rzu <- A - B1 * B2
  Su_tilde <- B1 / (sqrt(K) * Rxz * Rxsu - K * Rzu)
  Su_tilde * Sy
}


positive_beta <- function(Sigma, Rxsu, K){
  Sz <- sqrt(Sigma[3,3])
  Sxz <- Sigma[1,3]
  Rzu_facet <- facet_val(outer(Rxsu, K,
                               function(Rxsu, K) get_Rzu(Sigma, Rxsu, K)))
  Su_facet <- facet_val(outer(Rxsu, K,
                              function(Rxsu, K) get_Su(Sigma, Rxsu, K)))
  beta_facet <- getBetaIV(Sigma) - Su_facet * Rzu_facet * Sz / Sxz
  beta_facet > 0
}

in_prior <- function(Sigma, Rxsu, K, prior){
  Rzu <- facet_val(outer(Rxsu, K,
                         function(Rxsu, K) get_Rzu(Sigma, Rxsu, K)))
  Rzu_in_prior <- (Rzu > min(prior$Rzu)) & (Rzu < max(prior$Rzu))
  K <- facet_val(outer(rep(1, length(K)), K))
  K_in_prior <- (K > min(prior$K)) & (K < max(prior$K))
  Rxsu <- facet_val(outer(Rxsu, rep(1, length(Rxsu))))
  Rxsu_in_prior <- (Rxsu > min(prior$Rxsu)) & (Rxsu < max(prior$Rxsu))
  return(Rzu_in_prior & Rxsu_in_prior & K_in_prior)
}
