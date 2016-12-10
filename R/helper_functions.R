#' Construct vectors of points that outline a rectangle.
#'
#' @param xlower
#' @param xupper
#' @param ylower
#' @param yupper
#' @param step_x
#' @param step_y
#'
#' @return
#' @export
#'
#' @examples
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
#'
#' @return Matrix
#' @export
#'
#' @examples
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


#' Find the roots of a * x^2 + b * x + c = 0
#'
#' @param a Numeric vector
#' @param b Numeric vector
#' @param c Numeric vector
#'
#' @return Dataframe with columns x1, x2 and real. The first two columns contain
#' the complex-valued roots. The third column indicates if the roots are real.
#' @export
#'
#' @details This function computes an explicit solution unlike the base R
#' function \code{polyroot} which uses an iterative approach. Unlike
#' \code{polyroot}, \code{solve_quadratic} is vectorized in its arguments, but
#' only accepts real-valued coefficients as inputs.
#'
#' @examples
solve_quadratic <- function(a, b, c) {
  stopifnot(all(a != 0))
  stopifnot(is.atomic(a) && is.atomic(b) && is.atomic(c))
  stopifnot(is.numeric(a) && is.numeric(b) && is.numeric(c))
  n <- length(a)
  stopifnot(identical(n, length(b)))
  stopifnot(identical(n, length(c)))
  d <- b ^ 2 - 4 * a * c
  q <- rep(NA_complex_, n)
  real <- d >= 0
  # sign(0) = 0, but here we want it to be 1
  b_sign <- -1 * (b < 0) + 1 * (b >= 0)
  q[real] <- -0.5 * (b[real] + b_sign[real] * sqrt(d[real]))
  q[!real] <- -0.5 * (b[!real] + b_sign[!real] * 1i * sqrt(abs(d[!real])))
  x1 <- q / a # Function has already terminated if any(a == 0)
  ok <- q != 0 # q equals zero when b = c = 0
  x2 <- rep(NA_complex_, n)
  x2[ok] <- c[ok] / q[ok]
  x2[!ok] <- 0
  data.frame(root1 = x1, root2 = x2, real = real)
}

#' Find the real roots of a * x^2 + b * x + c = 0
#'
#' @param a Numeric Vector.
#' @param b Numeric Vector.
#' @param c Numeric Vector.
#'
#' @return Dataframe with two columns, containing the roots if real and NA
#' otherwise.
#' @export
#'
#' @details \code{solve_quadratic_real} simply calls \code{solve_quadratic} and
#' processes the results to remove any complex-valued roots. As such, it too is
#' vectorized in its arguments.
#'
#' @examples
solve_quadratic_real <- function(a, b, c){
  solution <- solve_quadratic(a, b, c)
  ok <- solution$real
  x1 <- x2 <- rep(NA_real_, length(ok))
  x1[ok] <- as.numeric(solution$root1[ok])
  x2[ok] <- as.numeric(solution$root2[ok])
  data.frame(x1, x2)
}

#' Find the roots of x^3 + a * x^2 + b * x + c = 0
#'
#' @param a Numeric vector.
#' @param b Numeric vector.
#' @param c Numeric vector.
#'
#' @return Dataframe with four columns: \code{x1}, \code{x2}, \code{x3},
#' and \code{three_real}. The column \code{x1} is real-valued since one root of
#' a cubic is necessarily real. The columns \code{x2} and \code{x3} are
#' complex-valued. Finally, \code{three_real} indicates whether the three roots
#' are all real-valued, i.e. if the imaginary parts of \code{x2} and \code{x3}
#' are both zero.
#' @details \code{solve_cubic_1} is not intended to be called directly: it is
#' used by wrapper function \code{solve_cubic} which calculates the roots of an
#' arbitrary cubic equation, i.e. one with no restriction on the first
#' coefficient. Unlike \code{polyroot}, \code{solve_cubic}only accepts
#' real-valued coefficients, is fully vectorized, and uses an explicit rather
#' than iterative solution method.
#' @export
#'
#' @examples
solve_cubic_1 <- function(a, b, c){
  stopifnot(is.atomic(a) && is.atomic(b) && is.atomic(c))
  stopifnot(is.numeric(a) && is.numeric(b) && is.numeric(c))
  n <- length(a)
  stopifnot(identical(n, length(b)))
  stopifnot(identical(n, length(c)))
  Q <- (a ^ 2 - 3 * b) / 9
  R <- (2 * a ^ 3 - 9 * a * b + 27 * c) / 54
  r3 <- R ^ 2 < Q ^ 3
  theta <- acos(R[r3] / sqrt(Q[r3] ^ 3))
  x1 <- rep(NA_real_, n)
  x2 <- x3 <- rep(NA_complex_, n)
  x1[r3] <- -2 * sqrt(Q[r3]) * cos(theta / 3) - a[r3] / 3
  x2[r3] <- -2 * sqrt(Q[r3]) * cos((theta + 2 * pi) / 3) - a[r3] / 3
  x3[r3] <- -2 * sqrt(Q[r3]) * cos((theta - 2 * pi) / 3) - a[r3] / 3
  A <- -1 * sign(R[!r3]) * (abs(R[!r3]) + sqrt(R[!r3] ^ 2 - Q[!r3] ^ 3)) ^ (1 / 3)
  ok <- A != 0
  B <- rep(NA_real_, length(A))
  B[ok] <- Q[!r3][ok] / A[ok]
  B[!ok] <- 0
  x1[!r3] <- (A + B) - a[!r3] / 3
  x2[!r3] <- -0.5 * (A + B) - (a[!r3] / 3) + 1i * (sqrt(3) / 2) * (A - B)
  x3[!r3] <- -0.5 * (A + B) - (a[!r3] / 3) - 1i * (sqrt(3) / 2) * (A - B)
  r3[!r3][abs(A - B) < sqrt(.Machine$double.eps)] <- TRUE
  data.frame(root1 = x1, root2 = x2, root3 = x3, three_real = r3)
}



#' Find the roots of a * x^3 + b * x^2 + c * x + d = 0
#'
#' @param a Numeric vector.
#' @param b Numeric vector.
#' @param c Numeric vector.
#' @param d Numeric vector.
#'
#' @return Dataframe with four columns: \code{x1}, \code{x2}, \code{x3}, and
#' \code{three_real}. The first three contain the roots while the third
#' indicates whether all three roots are real, as described in the help file for
#' \code{solve_cubic_1}.
#' @details This function is a wrapper to \code{solve_cubic_1} that allows the
#' coefficient on the cubic term to take on an arbitrary value.
#' @export
#'
#' @examples
solve_cubic <- function(a, b, c, d){
  stopifnot(all(a != 0))
  solve_cubic_1(b / a, c / a, d / a)
}

#' Find the real roots of a * x^3 + b * x^2 + c * x + d = 0
#'
#' @param a Numeric vector.
#' @param b Numeric vector.
#' @param c Numeric vector.
#' @param d Numeric vector.
#'
#' @return Dataframe containing the roots if real, NA otherwise.
#' @export
#'
#' @examples
solve_cubic_real <- function(a, b, c, d){
  solution <- solve_cubic(a, b, c, d)
  x1 <- solution$root1
  x2 <- x3 <- rep(NA_real_, length(x1))
  ok <- solution$three_real
  x2[ok] <- as.numeric(solution$root2[ok])
  x3[ok] <- as.numeric(solution$root3[ok])
  data.frame(root1 = x1, root2 = x2, root3 = x3)
}
