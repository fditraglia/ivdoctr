context("Example Datasets")

test_that("Replicate results of Acemoglu et al (2001)", {
  OLS1 <- lm(logpgp95 ~ avexpr, colonial)
  b_OLS1 <- as.numeric(round(coefficients(OLS1), 2)[2])
  expect_equal(0.52, b_OLS1)
  #OLS2 <- lm(logpgp95 ~ avexpr + lat_abst, colonial)
  #IV1 <- sem::tsls(logpgp95 ~ avexpr, ~ logem4, colonial)
  #IV2 <- sem::tsls(logpgp95 ~ avexpr + lat_abst, ~ logem4 + lat_abst, colonial)
})
