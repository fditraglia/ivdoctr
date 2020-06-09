#' @useDynLib ivdoctr, .registration = TRUE
NULL

#' Construct vectors of points that outline a rectangle.
#'
#' @param xleft The left side of the rectangle
#' @param xright The right side of the rectangle
#' @param ybottom The bottom of the rectangle
#' @param ytop The top of the rectangle
#' @param step_x The step size of the x coordinates
#' @param step_y The step size of the y coordinates
#'
#' @return List of x-coordinates and y-coordinates tracing the points around
#'   the rectangle
rect_points <- function(xleft, ybottom, xright, ytop, step_x, step_y){
  x_seq <- seq(xleft, xright, step_x)
  y_seq <- seq(ybottom, ytop, step_y)
  n_x <- length(x_seq)
  n_y <- length(y_seq)
  x <- c(x_seq, rep(x_seq[n_x], n_y - 2), rev(x_seq), rep(x_seq[1], n_y - 1))
  y <- c(rep(y_seq[1], n_x - 1), y_seq, rep(y_seq[n_y], n_x - 2), rev(y_seq))
  list(x = x, y = y)
}

#' Collapse 3-d array to matrix
#'
#' @param myarray A three-dimensional array.
#' @return Matrix with the 3rd dimension appended as rows to the matrix
collapse_3d_array <- function(myarray){
  out <- aperm(myarray, c(1, 3, 2))
  dim(out) <- c(dim(myarray)[1] * dim(myarray)[3], dim(myarray)[2])
  return(out)
}

#' Convert 3-d array to list of matrixes
#'
#' @param myArray A three-dimensional numeric array.
#' @return A list of numeric matrices.
#' @examples
#' M <- array(c(1, 1, 1, 1, 2, 2, 2, 2), c(2, 2, 2))
#' toList(M)
toList <- function(myArray){
  lapply(seq_len(dim(myArray)[3]), function(i) myArray[, , i, drop = TRUE])
}

#' Compute the share of draws that could contain a valid instrument.
#'
#' @param draws List of simulated draws
#' @return Numeric of the share of valid draws as determined by having the
#'   the restricted bounds for r_uz contain zero.
get_p_valid <- function(draws) {
  ans <- sum(draws$empty == FALSE & draws$restricted$r_uz_lower <= 0 &
             draws$restricted$r_uz_upper >= 0) / length(draws$empty)
  return(ans)
}
