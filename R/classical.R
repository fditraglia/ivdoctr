samplePosteriorClassical <- function(y_name, x_name, z_name,
                                     controls = NULL, data,
                                     drawSigma = covCLT, prior = NULL,
                                     n_Sigma_draws = 1e3,
                                     n_M_max_draws = 1e4,
                                     n_IdentSet_draws = 1e3, weight_sampling = TRUE){

  #Check Inputs
  stopifnot(names(formals(drawSigma)) == c("inData", "n_draws"))

  if(!is.null(prior)){
    stopifnot(all(names(prior) %in% c("K", "Rxsu", "Rzu")))
    stopifnot(length(prior$K) == 2)
    stopifnot(all((prior$K > 0) && (prior$K <= 1)))
    stopifnot(length(prior$Rxsu) == 2)
    stopifnot(all((prior$Rxsu > -1) && (prior$Rxsu < 1)))
    stopifnot(length(prior$Rzu) == 2)
    stopifnot(all((prior$Rzu > -1) && (prior$Rzu < 1)))
  }

  #Project out control regressors if present
  if(!is.null(controls)){
    y <- lm(reformulate(controls, response = y_name), data)$residuals
    xreg <- lm(reformulate(controls, response = x_name), data)
    x <- xreg$residuals
    xRsq <- summary(xreg)$r.squared
    zreg <- lm(reformulate(controls, response = z_name), data)
    z <- zreg$residuals
    zRsq <- summary(zreg)$r.squared
  }else{
    y <- get(y_name, data)
    x <- get(x_name, data)
    z <- get(z_name, data)
    xRsq <- zRsq <- 0 #No controls is the same as controls that are
                      #uncorrelated with the regressor and instrument
  }
  Sigma_draws <- drawSigma(inData = data.frame(x=x, y=y, z=z), n_Sigma_draws)
  Sigma_draws <- toList(Sigma_draws)
  Rho_draws <- lapply(Sigma_draws, cov2cor)
  Rho_draws_vech <- do.call(cbind, lapply(Rho_draws, vech))
  Sigma_MLE <- cov(cbind(x,y,z))
  Rho_MLE <- cov2cor(Sigma_MLE)

  if(is.null(prior)){
    K_L <- 0.01
    K_U <- 1
    Rxsu_L <- -0.99
    Rxsu_U <- 0.99
    Rzu_L <- -0.99
    Rzu_U <- 0.99
  }else{
    Rxsu_L <- prior$Rxsu[1]
    Rxsu_U <- prior$Rxsu[2]
    Rzu_L <- prior$Rzu[1]
    Rzu_U <- prior$Rzu[2]
    K_L <- prior$K[1]
    K_U <- prior$K[2]
  }

  # Convert prior bounds on Kappa to prior bounds on Kappa_tilde
  Ktilde_L <- max(xRsq + 0.01, (K_L - xRsq) / (1 - xRsq))
  Ktilde_U <- max(1, (K_U - xRsq) / (1 - xRsq))
  # Convert prior bounds on Rzu to prior bounds on Rzu_tilde
  RzuTilde_L <- Rzu_L / sqrt(1 - zRsq)
  RzuTilde_U <- Rzu_U / sqrt(1 - zRsq)
  if((RzuTilde_L < -1) | (RzuTilde_L > 1)) RzuTilde_L <- -0.99
  if((RzuTilde_U < -1) | (RzuTilde_U > 1)) RzuTilde_U <- 0.99
  # Convert prior bounds on Rxsu to prior bounds on Rxsu_tilde
  RxsuTilde_L <- Rxsu_L / sqrt(1 - xRsq / (Ktilde_U * (1 - xRsq) + xRsq))
  RxsuTilde_U <- Rxsu_U / sqrt(1 - xRsq / (Ktilde_L * (1 - xRsq) + xRsq))
  if((RxsuTilde_L < -1) | (RxsuTilde_L > 1)) RxsuTilde_L <- -0.99
  if((RxsuTilde_U < -1) | (RxsuTilde_U > 1)) RxsuTilde_U <- 0.99


  MLE <- classicalSampler(vech(Rho_MLE), n_IdentSet_draws, n_M_max_draws,
                          K_L, K_U, Rxsu_L, Rxsu_U, Rzu_L, Rzu_U, xRsq,weight_sampling)

  sim <- classicalSampler(Rho_draws_vech, n_IdentSet_draws, n_M_max_draws,
                          K_L, K_U, Rxsu_L, Rxsu_U, Rzu_L, Rzu_U, xRsq,weight_sampling)


  #--------------------------------------------------------------
  # Diagnose problems with first step (rejection sampler) at MLE
  #--------------------------------------------------------------
  # A value of 0 for step1eff indicates that the rejection sampler terminated
  #  by reaching the maximum allowable number of iterations before accepting
  #  the desired number of draws.
  if(MLE$step1eff == 0) stop("Rejection sampler efficiency < 10% at MLE")


  # Format MLE output more nicely
  # and add additional quantities of interest
  MLE$Sigma <- Sigma_MLE
  MLE$draws <- with(MLE, data.frame(K, Rxsu, Rzu, SuTilde, Ruv))
  MLE$K <- MLE$Rxsu <- MLE$Rzu <- MLE$SuTilde <- MLE$Ruv <- NULL
  MLE$draws$Su <- with(MLE, SuTilde2Su(Sigma, draws$SuTilde))
  MLE$draws$Beta <- with(MLE, getBeta(Sigma, draws$Su, draws$Rzu))
  MLE$draws$SuTilde <- MLE$draws$Ruv <- NULL
  MLE$diagnostic <- with(MLE, data.frame(Klower, Rxsu_Upper, step1eff, maxM))
  MLE$Klower <- MLE$step1eff <- MLE$maxM <- MLE$Rxsu_Upper <- NULL
  MLE$xRsq <- xRsq
  MLE$zRsq <- zRsq
  # Convert RxsuTilde, RzuTilde, Ktilde to non-tilde versions
  MLE$draws$K <- toKappa(MLE$draws$K, MLE$xRsq)
  MLE$draws$Rzu <- toRzu(MLE$draws$Rzu, MLE$zRsq)
  MLE$draws$Rxsu <- toRxsu(MLE$draws$Rxsu, MLE$xRsq, MLE$draws$K) #Not Ktilde!
  # Zeros to NAs if identified set is empty
  if(MLE$diagnostic$maxM == 0) MLE$draws[] <- NA

  # Format sim output more nicely
  # and add additional quantities of interest
  sim_diagnostic <- with(sim, data.frame(Klower, Rxsu_Upper, step1eff, maxM))
  sim$Klower <- sim$step1eff <- sim$maxM <- sim$Rxsu_Upper <- NULL
  sim <- lapply(sim, toList)
  sim$Sigma <- Sigma_draws
  sim$Su <- with(sim, Map(SuTilde2Su, Sigma, SuTilde))
  sim$Beta <- with(sim, Map(getBeta, Sigma, Su, Rzu))
  # Convert RxsuTilde, RzuTilde, Ktilde to non-tilde versions
  sim$K <- with(sim, Map(toKappa, K, xRsq))
  sim$Rzu <- with(sim, Map(toRzu, Rzu, zRsq))
  sim$Rxsu <- with(sim, Map(toRxsu, Rxsu, xRsq, K)) #Not Ktilde!
  # Zeros to NAs if identified set is empty
  empty_set <- which(sim_diagnostic$maxM == 0)
  sim$Sigma <- NULL #don't want to overwrite these with NAs! Remove, re-attach
  if(length(empty_set) > 0){
    for(i in 1:length(sim)){
      for(j in 1:length(empty_set)){
        sim[[i]][[empty_set[j]]] <- rep(NA, n_IdentSet_draws)
      }
    }
  }
  sim$Sigma <- Sigma_draws
  sim$draws <- with(sim, data.frame(K = unlist(K), Rxsu = unlist(Rxsu),
                                    Rzu = unlist(Rzu), Su = unlist(Su),
                                    Beta = unlist(Beta)))
  sim$draws_list <- with(sim, list(K = K, Rxsu = Rxsu, Rzu = Rzu,
                                   Su = Su, Beta = Beta))
  sim$diagnostic <- sim_diagnostic
  sim$Rzu <- sim$Rxsu <- sim$K <- sim$SuTilde <-
    sim$Ruv <- sim$Su <- sim$Beta <- NULL





  #-----------------------------------------------------------------------
  # Diagnose problems with first step (rejection sampler) at Sigma_draws
  #-----------------------------------------------------------------------
  # A value of 0 for step1eff indicates that the rejection sampler terminated
  #  by reaching the maximum allowable number of iterations before accepting
  #  the desired number of draws.
  step1fail_sim <- which(sim$step1eff == 0)
  if(length(step1fail_sim) > 0){
    warning(paste("Rejection sampler efficiency < 10% at",
                  length(step1fail_sim), "draws for Sigma"))
    fail_Sigma <- sim$Sigma[step1fail_sim]
    fail_diagnostic <- sim$diagnostic[step1fail_sim,]
    sim$fail <- list(Sigma = fail_Sigma, diagnostic = fail_diagnostic)
  }

  return(list(MLE = MLE, sim = sim))

}

#============================== Functions for plotting the results

plot.full.classical <- function(Sigma, xRsq, zRsq, prior = NULL, theta, phi,
                                TeX = FALSE){
  R <- cov2cor(Sigma)
  Rxy <- R[1,2]
  Rxz <- R[1,3]
  Rzy <- R[2,3]
  Rxsu <- seq(-0.99, 0.99, length.out = 50)
  K <- seq(((Rxy^2)+(Rxz^2)-2*Rxy*Rxz*Rzy)/(1-(Rzy^2)), 1, length.out = 50)


  # Calculate RzuTilde, NA if it violates unit circle restriction
  # or it if lies outside (-1, 1)
  Rzu <- get_Rzu_matrix(Sigma, Rxsu, K)

  # We have been working with the transformed versions, i.e. "tilde"
  KTilde <- K
  RzuTilde <- Rzu
  RxsuTilde <- Rxsu

  # Transform Rzu, Rxsu, K to "non-tilde" versions
  K <- toKappa(KTilde, xRsq)
  Rzu <- toRzu(RzuTilde, zRsq)
  Rxsu <- toRxsu(RxsuTilde, xRsq, K) #Kappa without tilde!

  # Clean up!
  rm(KTilde, RzuTilde, RxsuTilde)

  if(is.null(prior)){
    colors <- "lightgreen"
  }else{
    colors <- ifelse(in_prior(Sigma, Rxsu, K, prior), "lightgreen", "indianred2")
  }


  if(TeX){
    x_lab <- "$\\rho_{T^*u}$"
    y_lab <- "$\\kappa$"
    z_lab <- "$\\rho_{uz}$"
  }else{
    x_lab <- "Cor(T*,u)"
    y_lab <- "Kappa"
    z_lab <- "Cor(z,u)"
  }

  persp(sort(Rxsu), sort(K), Rzu[order(Rxsu),order(K)],
        theta = theta, phi = phi, xlab = x_lab, ylab = y_lab,
        zlab = z_lab, ticktype = "detailed", col = colors, shade = 0.3,
        border = NA, bg = "white")
}


plot.pos.classical <- function(Sigma, xRsq, zRsq, prior = NULL, theta, phi,
                               TeX = FALSE){
  R <- cov2cor(Sigma)
  Rxy <- R[1,2]
  Rxz <- R[1,3]
  Rzy <- R[2,3]

  if(is.null(prior)){
    RxsuTilde <- seq(-0.99, 0.99, length.out = 50)
    KTilde <- seq(((Rxy^2)+(Rxz^2)-2*Rxy*Rxz*Rzy)/(1-(Rzy^2)), 1, length.out = 50)

  }else{
    prior$KTilde <- toKappaTilde(prior$K, xRsq)
    KTilde_max <- min(1, max(prior$KTilde))
    KTilde_min <- max(((Rxy^2)+(Rxz^2)-2*Rxy*Rxz*Rzy)/(1-(Rzy^2)), min(prior$KTilde))
    KTilde <- seq(KTilde_min, KTilde_max, length.out = 50)


    RxsuTilde_L <- toRxsuTilde(min(prior$Rxsu), xRsq, toKappa(KTilde_max, xRsq))
    RxsuTilde_U <- toRxsuTilde(max(prior$Rxsu), xRsq, toKappa(KTilde_min, xRsq))
    if((RxsuTilde_L < -1) | (RxsuTilde_L > 1)) RxsuTilde_L <- -0.99
    if((RxsuTilde_U < -1) | (RxsuTilde_U > 1)) RxsuTilde_U <- 0.99
    RxsuTilde <- seq(RxsuTilde_L, RxsuTilde_U, length.out = 50)
  }

  # Calculate RzuTilde, NA if it violates unit circle restriction
  # or it if lies outside (-1, 1)
  RzuTilde <- get_Rzu_matrix(Sigma, RxsuTilde, KTilde)

  # Color in region that maps to positive beta
  if(is.null(prior)){
    colors <- ifelse(positive_beta(Sigma, RxsuTilde, KTilde),
                     "dodgerblue2", "indianred2")
  }else{
    colors <- ifelse(positive_beta(Sigma, RxsuTilde, KTilde),
                     "dodgerblue2", "lightgreen")
  }

  if(TeX){
    x_lab <- "$\\rho_{T^*u}$"
    y_lab <- "$\\kappa$"
    z_lab <- "$\\rho_{uz}$"
  }else{
    x_lab <- "Cor(T*,u)"
    y_lab <- "Kappa"
    z_lab <- "Cor(z,u)"
  }

  # Convert to non-tilde
  Rzu <- toRzu(RzuTilde, zRsq)
  K <- toKappa(KTilde, xRsq)
  Rxsu <- toRxsu(RxsuTilde, xRsq, K) #Not KTilde!
  rm(RzuTilde, RxsuTilde, KTilde)

  # Plot non-tilde stuff: notice that we ensure ascending values for "x" and "y"
  # and arrange "z" accordingly. This makes sure persp doesn't get confused.
  # Note that our transformation for Rxsu can potentially give a non-monotonic
  # order although this is unlikely
  persp(sort(Rxsu), sort(K), Rzu[order(Rxsu),order(K)],
        theta = theta, phi = phi, xlab = x_lab, ylab = y_lab,
        zlab = z_lab, ticktype = "detailed", col = colors, shade = 0.3,
        border = NA, bg = "white")
}
