context("Datasets")

test_that("Replicate Table 4 of Acemoglu et al (2001)", {
  attach(colonial)
  b <- function(reg) as.numeric(round(coefficients(reg), 2)[2])
  SE <- function(reg) round(summary(reg)$coefficients[2,2], 2)
  OLS1 <- lm(logpgp95 ~ avexpr)
  OLS2 <- lm(logpgp95 ~ avexpr + lat_abst)
  IV1 <- AER::ivreg(logpgp95 ~ avexpr | logem4)
  IV2 <- AER::ivreg(logpgp95 ~ avexpr + lat_abst | logem4 + lat_abst)
  detach(colonial)
  expect_equal(0.52, b(OLS1))
  expect_equal(0.06, SE(OLS1))
  expect_equal(0.47, b(OLS2))
  expect_equal(0.06, SE(OLS2))
  expect_equal(0.94, b(IV1))
  expect_equal(0.16, SE(IV1))
  expect_equal(1, b(IV2))
  expect_equal(0.22, SE(IV2))
})

test_that("Replicate Tables II and III of Becker & Woessman (2009)", {
  attach(weber)
  b <- function(reg) as.numeric(round(coefficients(reg), 3)[2])
  SE <- function(reg) round(summary(reg)$coefficients[2,2], 3)
  controls <- c("f_young", "f_jew", "f_fem", "f_ortsgeb",
                "f_pruss", "hhsize", "lnpop", "gpop",
                "f_miss", "f_blind", "f_deaf", "f_dumb")
  second_stage <- reformulate(c("f_prot", controls), "f_rw")
  first_stage <- reformulate(c("kmwittenberg", controls))
  OLS1 <- lm(f_rw ~ f_prot)
  expect_equal(0.080, b(OLS1))
  expect_equal(0.015, SE(OLS1))
  OLS2 <- lm(second_stage)
  expect_equal(0.099, b(OLS2))
  expect_equal(0.010, SE(OLS2))
  IV <- AER::ivreg(f_rw ~ f_prot + f_young + f_jew + f_fem + f_ortsgeb + f_pruss +
                     hhsize + lnpop + gpop + f_miss + f_blind + f_deaf + f_dumb |
                     . - f_prot + kmwittenberg)
  expect_equal(0.189, b(IV))
  expect_equal(0.028, SE(IV))
  detach(weber)
})
