context("Helper Functions")

test_that("toList", {
  m1 <- rep(1, 4)
  m2 <- rep(2, 4)
  arrayM <- array(c(m1, m2), c(2, 2, 2))
  listM <- list(matrix(m1, 2, 2), matrix(m2, 2, 2))
  expect_equal(listM, toList(arrayM))
})

test_that("rinvwish agrees with riwish from MCMCpack", {
  skip_if_not_installed('MCMCpack')
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

test_that("get_p_valid correctly computes share of potentially valid instruments
          when bounds cover 0", {
  draws <- list(empty = c(1, 0, 0, 0),
                restricted = list(r_uz_lower = c(-1, -1, -1, -1),
                                  r_uz_upper = c(1, 1, 1, 1)))
  ans <- get_p_valid(draws)
  expected <- 0.75
})

test_that("get_p_valid correctly computes share of potentially valid instruments
          when some bounds cover zero", {
  draws <- list(empty = c(1, 0, 0, 0),
                restricted = list(r_uz_lower = c(-1, 0.5, -1, -1),
                                  r_uz_upper = c(1, 1, 1, 1)))
  ans <- get_p_valid(draws)
  expected <- 0.5
})

test_that("get_p_valid correctly computes share of potentially valid instruments
          when bounds cover 0, but all identified sets are empty", {
  draws <- list(empty = c(1, 1, 1, 1),
                restricted = list(r_uz_lower = c(-1, -1, -1, -1),
                                  r_uz_upper = c(1, 1, 1, 1)))
  ans <- get_p_valid(draws)
  expected <- 0
})

test_that("get_p_valid correctly computes share of potentially valid instruments
          when no bounds cover 0", {
  draws <- list(empty = c(0, 0, 0, 0),
                restricted = list(r_uz_lower = c(0.5, 0.5, 0.5, 0.5),
                                  r_uz_upper = c(1, 1, 1, 1)))
  ans <- get_p_valid(draws)
  expected <- 0
})

test_that("get_p_valid correctly computes share of potentially valid instruments
          when no bounds cover 0", {
  draws <- list(empty = c(0, 0, 0, 0),
                restricted = list(r_uz_lower = c(-1, 0.5, 0.5, -1),
                                  r_uz_upper = c(1, 1, 1, 1)))
  ans <- get_p_valid(draws)
  expected <- 0.5
})


