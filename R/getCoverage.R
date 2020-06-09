#' Computes coverage of list of intervals
#'
#' @param data 2-column data frame of confidence intervals
#' @param guess 2-element vector of confidence interval
#'
#' @return Coverage percentage
#'
getCoverage <- function(data, guess) {
  coverage <- (min(guess) <= apply(data, 1, min)) * (max(guess) >= apply(data, 1, max))
  coverage <- sum(coverage) / nrow(data)
  return(coverage)
}

#' Generates smallest covering interval
#'
#' @param data 2-column data frame of confidence intervals
#' @param center 2-element vector to center coverage interval
#' @param conf Confidence level
#' @param tol Tolerance level for convergence
#'
#' @return 2-element vector of confidence interval
#'
getInterval <- function(data, center, conf = 0.9, tol = 1e-6) {
  data <- na.omit(data)
  minInt <- center
  stretch <- max(abs(min(minInt) - min(data)),
                 abs(max(minInt) - max(data)))
  maxInt <- c(min(center) - stretch, max(center) + stretch)
  delta <- 1
  prev <- maxInt
  while (delta > tol) {
    eps <- abs(min(maxInt) - min(minInt)) / 2
    miss <- getCoverage(data, prev) - conf

    if (miss >= 0) {
      guess <- c(min(prev) + eps, max(prev) - eps)
    } else {
      guess <- c(min(prev - eps), max(prev) + eps)
    }

    coverage <- getCoverage(data, guess)

    if (coverage - conf >= 0) {
      delta <- abs(maxInt - guess)[1]
      maxInt <- guess
      prev <- guess
    } else {
      delta <- abs(minInt - guess)[1]
      minInt <- guess
      prev <- guess
    }
  }
  return(maxInt)
}
