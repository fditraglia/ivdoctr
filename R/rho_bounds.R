obs <- foo
r_TstarU <- 0.3

c0 <- with(obs, -0.25 * r_Ty^6 * r_Tz^2)

c1_1 <- with(obs, r_Ty^5 * r_Tz * r_zy)
c1_2 <- with(obs, r_Ty^4 * r_Tz^2)
c1 <- -0.5 * c1_1 + c1_2 + 0.25 * (2 * c1_1 - c1_2) * r_TstarU

c2_1 <- with(obs, r_Ty^4 * r_zy^2)
c2_2 <- with(obs, r_Ty^3 * r_Tz * r_zy)
c2_3 <- with(obs, r_Ty^2 * r_Tz^2)
c2 <- -0.25 * c2_1 + c2_2 - c2_3 + 0.25 * (c2_1 - 4 * c2_2 + c2_3) * r_TstarU

c3 <- with(obs, 0.25 * r_TstarU^2 * r_Tz^2)

coeffs <- c(c0, c1, c2, c3)
roots <- polyroot(coeffs)

solve_quadratic <- function(a, b, c) {
  stopifnot(is.atomic(a) && is.atomic(b) && is.atomic(c))
  stopifnot(is.numeric(a) && is.numeric(b) && is.numeric(c))
  n <- length(a)
  stopifnot(identical(n, length(b)))
  stopifnot(identical(n, length(c)))
  d <- b^2 - 4 * a * c
  q <- rep(NA_real_, length(a))
  real <- d >= 0
  q[real] <- -0.5 * (b[real] + sign(b[real]) * sqrt(d[real]))
  q[!real] <- -0.5 * (b[!real] + sign(b[!real]) * 1i * sqrt(abs(d[!real])))
  x1 <- q / a
  x2 <- c / q
  list(root1 = x1, root2 = x2)
}

solve_cubic_1 <- function(a, b, c){
  stopifnot(is.atomic(a) && is.atomic(b) && is.atomic(c))
  stopifnot(is.numeric(a) && is.numeric(b) && is.numeric(c))
  n <- length(a)
  stopifnot(identical(n, length(b)))
  stopifnot(identical(n, length(c)))
  n <- length(a)
  stopifnot(identical(n, length(b)))
  stopifnot(identical(n, length(c)))
  Q <- (a^2 - 3 * b) / 9
  R <- (2 * a^3 - 9 * a * b + 27 * c) / 54
  r3 <- R^2 < Q^3
  theta <- acos(R[all_real] / sqrt(Q[all_real]^3))
  x1[r3] <- -2 * sqrt(Q[r3]) * cos(theta / 3) - a[all_real] / 3
  x2[r3] <- -2 * sqrt(Q[r3]) * cos((theta + 2 * pi) / 3) - a[r3] / 3
  x3[r3] <- -2 * sqrt(Q[r3]) * cos((theta - 2 * pi) / 3) - a[r3] / 3
  A <- -1 * sign(R[!r3]) * (abs(R[!r3]) + sqrt(R[!r3]^2 - Q[!r3]^3))^(1 / 3)
  # B <-
}
















