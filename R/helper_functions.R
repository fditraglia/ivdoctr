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
  x <- as.matrix(inData$x)
  y <- as.matrix(inData$y)
  z <- as.matrix(inData$z)
  e1 <- lm.fit(x, y)$residuals # OLS resids
  e2 <- lm.fit(z, y)$residuals # IV first-stage resids
  xe1 <- x * e1
  ze2 <- z * e2
  Omega <- cov(cbind(xe1, ze2)) / length(x)
  Sigma <- cov(inData)
  sims <- MASS::mvrnorm(n_draws, c(cov(x,y), cov(z,y)), Omega)
  upper_det <- var(x) * var(y) - sims[,1]^2
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
