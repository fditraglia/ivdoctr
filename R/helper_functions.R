toList <- function(myArray){
#-------------------------------------------------------------
# convert 3d array of matrices to list of matrices
#-------------------------------------------------------------
  lapply(seq_len(dim(myArray)[3]), function(i) myArray[,,i,drop = TRUE])
}

covJeffreys <- function(inData, n_draws){
  n <- nrow(inData)
  S <- (n - 1) * cov(inData)
  return(rinvwish(n_draws, n - 1, S))
}


covCLT <- function(inData, n_draws){
  stopifnot(setequal(names(inData), c("x", "y", "z")))
  x <- as.matrix(inData$x)
  y <- as.matrix(inData$y)
  z <- as.matrix(inData$z)
  e1 <- lm.fit(x, y)$residuals
  e2 <- lm.fit(z, y)$residuals
  xe1 <- x * e1
  ze2 <- z * e2
  Omega <- cov(cbind(xe1, ze2)) / length(x)
  Sigma <- cov(inData)
  sims <- MASS::mvrnorm(n_draws, c(cov(x,y), cov(z,y)), Omega)
  g <- function(sim_row){
    out <- Sigma
    out[1,2] <- out[2,1] <- sim_row[1]
    out[2,3] <- out[3,2] <- sim_row[2]
    return(out)
  }
  Sigma_draws <- apply(sims, 1, g)
  dim(Sigma_draws) <- c(3, 3, n_draws)
  pd <- apply(Sigma_draws, 3, function(M) all(eigen(M)$values > 0))
  stopifnot(all(pd))
  return(Sigma_draws)
}
