context("Helper Functions")

test_that("quadratic solvers",{
  a <- c(1, 2, 2, 4)
  b <- c(0, 0, 0, -16)
  c <- c(-1, 0, 1, 16)
  expected_answer1 <- data.frame(x1 = c(-1, 0, NA, 2), x2 = c(1, 0, NA, 2))
  ans1 <- solve_quadratic_real(a, b, c)
  expect_equal(expected_answer1, ans1)
  
  expected_answer2 <- data.frame(root1 = -1i/sqrt(2), root2 = 1i/sqrt(2), real = FALSE)
  ans2 <- solve_quadratic(2, 0, 1)
  expect_equal(expected_answer2, ans2)
})

test_that("cubic solvers", {
  expect_equal(solve_cubic(1, 0, 0, 0),
               data.frame(root1 = 0, root2 = 0 + 0i, root3 = 0 + 0i,
                          three_real = TRUE))
})

test_that("vech and devech", {
  M <- matrix(c(11, 12, 13, 14,
                12, 22, 23, 24,
                13, 23, 33, 34,
                14, 24, 34, 44), 4, 4, byrow = TRUE)
  v <- drop(vech(M))
  expect_equal(v, c(11:14, 22:24, 33:34, 44))
  expect_that(vech(matrix(1:6, 3, 2)), throws_error())
  expect_equal(devech(v, 4), M)
  expect_that(devech(v, 3), throws_error())
})

test_that("toList", {
  m1 <- rep(1, 4)
  m2 <- rep(2, 4)
  arrayM <- array(c(m1, m2), c(2, 2, 2))
  listM <- list(matrix(m1, 2, 2), matrix(m2, 2, 2))
  expect_equal(listM, toList(arrayM))
})

test_that("rinvwish agrees with riwish from MCMCpack", {
  V <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  set.seed(4321)
  riwish_draw1 <- MCMCpack::riwish(10, 10 * V)
  riwish_draw2 <- MCMCpack::riwish(10, 10 * V)
  riwish_draw <- array(c(riwish_draw1, riwish_draw2), c(2, 2, 2))
  set.seed(4321)
  rinvwish_draw <- drop(rinvwish(2, 10, 10 * V))
  expect_equal(riwish_draw, rinvwish_draw)
})

test_that("error on improper input to covCLT", {
  bad_input1 <- data.frame(x = rnorm(3), y = rnorm(3))
  bad_input2 <- data.frame(q = rnorm(3), y = rnorm(3))
  expect_that(covCLT(bad_input1, 10), throws_error())
  expect_that(covCLT(bad_input2, 10), throws_error())
})

test_that("error on non-pd cov matrix draw from covCLT", {
  set.seed(1234)
  x <- rnorm(10)
  y <- x + rnorm(10) / 100
  z <- x + rnorm(10) / 100
  foo <- data.frame(x, y, z)
  expect_that(covCLT(foo, 100), throws_error())
})

test_that("rectangle points are generated properly", {
  xleft <- 0
  ybottom <- 0
  xright <- 1
  ytop <- 1
  step_x <- 1
  step_y <- 1
  expected <- list(x = c(0, 1, 1, 0, 0),
                   y = c(0, 0, 1, 1, 0))
  ans <- rect_points(xleft, ybottom, xright, ytop, step_x, step_y)
  expect_equal(expected, ans)
})

test_that("rectangle points are generated properly", {
  xleft <- 0
  ybottom <- 0
  xright <- 1
  ytop <- 1
  step_x <- 0.5
  step_y <- 0.5
  expected <- list(x = c(0, 0.5, 1, 1, 1, 0.5, 0, 0, 0),
                   y = c(0, 0, 0, 0.5, 1, 1, 1, 0.5, 0))
  ans <- rect_points(xleft, ybottom, xright, ytop, step_x, step_y)
  expect_equal(expected, ans)
})
