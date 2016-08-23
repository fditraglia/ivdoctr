get_estimates <- function(y_name, T_name, z_name, data, controls = NULL,
                          robust = FALSE) {
  first_stage <- reformulate(c(z_name, controls), response = NULL)
  second_stage <- reformulate(c(T_name, controls), response = y_name)
  OLS <- lm(second_stage, data)
  IV <-  AER::ivreg(second_stage, first_stage, data)

  b_OLS <- coef(OLS)[T_name]
  b_IV <- coef(IV)[T_name]

  if (robust) {
    se_OLS <- sqrt(diag(sandwich::vcovHC(OLS, type = 'HC0')))[T_name]
    se_IV <- sqrt(diag(sandwich::vcovHC(IV, type = 'HC0')))[T_name]
  } else {
    se_OLS <- sqrt(diag(vcov(OLS)))[T_name]
    se_IV <- sqrt(diag(vcov(IV)))[T_name]
  }

  list(n = nrow(data),
       b_OLS = b_OLS,
       se_OLS = se_OLS,
       b_IV = b_IV,
       se_IV = se_IV)
}

get_HPDI <- function(draws, level = 0.9) {
  interval <- coda::HPDinterval(coda::as.mcmc(draws), level)
  lower <- interval[[1]]
  upper<- interval[[2]]
  return(data.frame(lower = lower, median = median(draws), upper = upper))
}

summarize_bounds <- function(draws) {
  unrestricted <- with(draws$unrestricted,
                       rbind(k_lower = get_HPDI(k_lower),
                             r_uz_lower = get_HPDI(r_uz_lower),
                             r_uz_upper = get_HPDI(r_uz_upper)))
  # Add in restricted r_uz_lower and upper bounds later:
  restricted <- with(draws$restricted,
                     rbind(beta_lower = get_HPDI(beta_lower),
                           beta_upper = get_HPDI(beta_upper)))
  list(unrestricted = unrestricted,
       r_TstarU_restriction = draws$r_TstarU_restriction,
       k_restriction = draws$k_restriction,
       p_empty = mean(draws$empty),
       restricted = restricted)
}

summarize_posterior <- function(draws) {
  HPDI <- with(draws$posterior, rbind(r_uz = get_HPDI(r_uz),
                                      beta = get_HPDI(beta)))
  list(r_TstarU_restriction = draws$r_TstarU_restriction,
       k_restriction = draws$k_restriction,
       p_empty = mean(draws$empty),
       HPDI = HPDI)
}


# #============================== Functions for plotting the results
#
# plot.full.classical <- function(Sigma, xRsq, zRsq, prior = NULL, theta, phi,
#                                 TeX = FALSE){
#   R <- cov2cor(Sigma)
#   Rxy <- R[1,2]
#   Rxz <- R[1,3]
#   Rzy <- R[2,3]
#   Rxsu <- seq(-0.99, 0.99, length.out = 50)
#   K <- seq(((Rxy^2)+(Rxz^2)-2*Rxy*Rxz*Rzy)/(1-(Rzy^2)), 1, length.out = 50)
#
#
#   # Calculate RzuTilde, NA if it violates unit circle restriction
#   # or it if lies outside (-1, 1)
#   Rzu <- get_Rzu_matrix(Sigma, Rxsu, K)
#
#   # We have been working with the transformed versions, i.e. "tilde"
#   KTilde <- K
#   RzuTilde <- Rzu
#   RxsuTilde <- Rxsu
#
#   # Transform Rzu, Rxsu, K to "non-tilde" versions
#   K <- toKappa(KTilde, xRsq)
#   Rzu <- RzuTilde ## toRzu(RzuTilde, zRsq)
#   Rxsu <- RxsuTilde ##toRxsu(RxsuTilde, xRsq, K) #Kappa without tilde!
#
#   # Clean up!
#   rm(KTilde, RzuTilde, RxsuTilde)
#
#   if(is.null(prior)){
#     colors <- "lightgreen"
#   }else{
#     colors <- ifelse(in_prior(Sigma, Rxsu, K, prior), "lightgreen", "indianred2")
#   }
#
#
#   if(TeX){
#     x_lab <- "$\\rho_{T^*u}$"
#     y_lab <- "$\\kappa$"
#     z_lab <- "$\\rho_{uz}$"
#   }else{
#     x_lab <- "Cor(T*,u)"
#     y_lab <- "Kappa"
#     z_lab <- "Cor(z,u)"
#   }
#
#   persp(sort(Rxsu), sort(K), Rzu[order(Rxsu),order(K)],
#         theta = theta, phi = phi, xlab = x_lab, ylab = y_lab,
#         zlab = z_lab, ticktype = "detailed", col = colors, shade = 0.3,
#         border = NA, bg = "white")
# }
#
#
# plot.pos.classical <- function(Sigma, xRsq, zRsq, prior = NULL, theta, phi,
#                                TeX = FALSE){
#   R <- cov2cor(Sigma)
#   Rxy <- R[1,2]
#   Rxz <- R[1,3]
#   Rzy <- R[2,3]
#
#   if(is.null(prior)){
#     RxsuTilde <- seq(-0.99, 0.99, length.out = 50)
#     KTilde <- seq(((Rxy^2)+(Rxz^2)-2*Rxy*Rxz*Rzy)/(1-(Rzy^2)), 1, length.out = 50)
#
#   }else{
#     prior$KTilde <- toKappaTilde(prior$K, xRsq)
#     KTilde_max <- min(1, max(prior$KTilde))
#     KTilde_min <- max(((Rxy^2)+(Rxz^2)-2*Rxy*Rxz*Rzy)/(1-(Rzy^2)), min(prior$KTilde))
#     KTilde <- seq(KTilde_min, KTilde_max, length.out = 50)
#
#
#     RxsuTilde_L <- min(prior$Rxsu) ## toRxsuTilde(min(prior$Rxsu), xRsq, toKappa(KTilde_max, xRsq))
#     RxsuTilde_U <- max(prior$Rxsu) ## toRxsuTilde(max(prior$Rxsu), xRsq, toKappa(KTilde_min, xRsq))
#     if((RxsuTilde_L < -1) | (RxsuTilde_L > 1)) RxsuTilde_L <- -0.99
#     if((RxsuTilde_U < -1) | (RxsuTilde_U > 1)) RxsuTilde_U <- 0.99
#     RxsuTilde <- seq(RxsuTilde_L, RxsuTilde_U, length.out = 50)
#   }
#
#   # Calculate RzuTilde, NA if it violates unit circle restriction
#   # or it if lies outside (-1, 1)
#   RzuTilde <- get_Rzu_matrix(Sigma, RxsuTilde, KTilde)
#
#   # Color in region that maps to positive beta
#   if(is.null(prior)){
#     colors <- ifelse(positive_beta(Sigma, RxsuTilde, KTilde),
#                      "dodgerblue2", "indianred2")
#   }else{
#     colors <- ifelse(positive_beta(Sigma, RxsuTilde, KTilde),
#                      "dodgerblue2", "lightgreen")
#   }
#
#   if(TeX){
#     x_lab <- "$\\rho_{T^*u}$"
#     y_lab <- "$\\kappa$"
#     z_lab <- "$\\rho_{uz}$"
#   }else{
#     x_lab <- "Cor(T*,u)"
#     y_lab <- "Kappa"
#     z_lab <- "Cor(z,u)"
#   }
#
#   # Convert to non-tilde
#   Rzu <- RzuTilde ## toRzu(RzuTilde, zRsq)
#   K <- toKappa(KTilde, xRsq)
#   Rxsu <- RxsuTilde ## toRxsu(RxsuTilde, xRsq, K) #Not KTilde!
#   rm(RzuTilde, RxsuTilde, KTilde)
#
#   # Plot non-tilde stuff: notice that we ensure ascending values for "x" and "y"
#   # and arrange "z" accordingly. This makes sure persp doesn't get confused.
#   # Note that our transformation for Rxsu can potentially give a non-monotonic
#   # order although this is unlikely
#   persp(sort(Rxsu), sort(K), Rzu[order(Rxsu),order(K)],
#         theta = theta, phi = phi, xlab = x_lab, ylab = y_lab,
#         zlab = z_lab, ticktype = "detailed", col = colors, shade = 0.3,
#         border = NA, bg = "white", 0.75, 1.25)
# }
