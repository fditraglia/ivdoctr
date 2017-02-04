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
