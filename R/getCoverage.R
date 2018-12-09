#' Computes coverage of list of intervals
#'
#' @param data 2-column data frame of confidence intervals
#' @param guess 2-element vector of confidence interval
#'
#' @return Coverage percentage
#' @export
#'
#' @examples
getCoverage <- function(data, guess) {
  coverage <- (min(guess) <= data$min) * (max(guess) >= data$max)
  coverage <- sum(coverage) / nrow(data)
  return(coverage)
}

#' Generates smallest covering interval
#'
#' @param data 2-column data frame of confidence intervals
#' @param conf Confidence level
#' @param tol Tolerance level for convergence
#'
#' @return 2-element vector of confidence interval
#' @export
#'
#' @examples
getInterval <- function(data, conf = 0.95, tol = 1e-6) {
  maxInt <- c(min(data), max(data))
  minInt <- c(mean(maxInt), mean(maxInt))
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
