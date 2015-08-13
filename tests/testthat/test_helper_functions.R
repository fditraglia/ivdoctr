context("Helper Functions")

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
