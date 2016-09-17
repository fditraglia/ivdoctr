context("Identified Set")

test_that("get_k_tilde_lower computes properly", {
          obs <- list(r_Ty = 0.5,
                      r_Tz = 0.5,
                      r_zy = 0.5)
          expect_equal(get_k_tilde_lower(obs), 1/3)
})

test_that("get_k_lower computes properly without controls", {
          obs <- list(r_Ty = 0.5,
                      r_Tz = 0.5,
                      r_zy = 0.5,
                      T_Rsq = 0)
          expect_equal(get_k_lower(obs), 1/3)
})

test_that("get_k_lower computes properly with controls", {
          obs <- list(r_Ty = 0.5,
                      r_Tz = 0.5,
                      r_zy = 0.5,
                      T_Rsq = 0.25)
          expect_equal(get_k_lower(obs), 1/2)
})
