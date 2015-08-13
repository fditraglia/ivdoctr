samplePosteriorClassical <- function(y_name, x_name, z_name,
                                     controls = NULL, data,
                                     drawSigma = covCLT, prior = NULL,
                                     n_Sigma_draws = 1e3,
                                     n_M_max_draws = 1e4,
                                     n_IdentSet_draws = 1e3){

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

  if(is.null(prior)){
    K_L <- 0.01
    K_U <- 1
    Rxsu_L <- -0.99
    Rxsu_U <- 0.99
    Rzu_L <- -0.99
    Rzu_U <- 0.99
  }else{
    K_L <- prior$K[1]
    K_U <- prior$K[2]
    Rxsu_L <- prior$Rxsu[1]
    Rxsu_U <- prior$Rxsu[2]
    Rzu_L <- prior$Rzu[1]
    Rzu_U <- prior$Rzu[2]
  }

  #Project out control regressors if present
  if(!is.null(controls)){
    y <- lm(reformulate(controls, response = y_name), data)$residuals
    xreg <- lm(reformulate(controls, response = x_name), data)
    x <- xreg$residuals
    xRsq <- summary(xreg)$r.squared
    z <- lm(reformulate(controls, response = z_name), data)$residuals
  }else{
    y <- get(y_name, data)
    x <- get(x_name, data)
    xRsq <- 0 # No controls is the same as controls that are uncorrelated
              # with the variable of interest
    z <- get(z_name, data)
  }
  Sigma_draws <- drawSigma(inData = data.frame(x=x, y=y, z=z), n_Sigma_draws)
  Sigma_draws <- toList(Sigma_draws)
  Rho_draws <- lapply(Sigma, cov2cor)
  Rho_draws_vech <- do.call(cbind, lapply(Rho, vech))
  Sigma_MLE <- cov(cbind(x,y,z))
  Rho_MLE <- cov2cor(Sigma_MLE)

  MLE <- classicalSampler(vech(Rho_MLE), n_IdentSet_draws, n_M_max_draws,
                          K_L, K_U, Rxsu_L, Rxsu_U, Rzu_L, Rzu_U)
  sim <- classicalSampler(Rho_draws_vech, n_IdentSet_draws, n_M_max_draws,
                          K_L, K_U, Rxsu_L, Rxsu_U, Rzu_L, Rzu_U)

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
  MLE$draws$Su <- with(MLE, draws$SuTilde / sqrt(Sigma[2,2]))
  MLE$draws$Szu <- with(MLE, draws$Rzu * draws$Su * sqrt(Sigma[3,3]))
  MLE$draws$Beta <- with(MLE, (Sigma[3,2] - draws$Su) / Sigma[3,1])
  MLE$draws$SuTilde <- MLE$draws$Szu <- MLE$draws$Ruv <- NULL
  MLE$diagnostic <- with(MLE, data.frame(Klower, step1eff, maxM))
  MLE$Klower <- MLE$step1eff <- MLE$maxM <- NULL


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

  }


}

getPosteriorSamples <- function(x, y, z,
                                n_Sigma_draws = 1e2,
                                n_IdentSet_draws = 1e3,
                                n_M_max_draws = 1e4,
                                K_bounds = c(0.2, 1),
                                Rxsu_bounds = c(-0.9, 0.9),
                                Rzu_bounds = c(-0.9, 0.9)){

  # If there are control variables, work with kappa_tilde
  # in the underlying C++ routine

  # add in some sanity checks / coercion on inputs
  inData <- cbind(x, y, z)
  posterior <- postProcess(samplePosterior(inData, n_Sigma_draws,
                              n_IdentSet_draws, n_M_max_draws,
                              K_bounds[1], K_bounds[2],
                              Rxsu_bounds[1], Rxsu_bounds[2],
                              Rzu_bounds[1], Rzu_bounds[2]))

  priorMLE <- postProcess(samplePosterior(inData, 1,
                              n_IdentSet_draws, n_M_max_draws,
                              K_bounds[1], K_bounds[2],
                              Rxsu_bounds[1], Rxsu_bounds[2],
                              Rzu_bounds[1], Rzu_bounds[2]))
  return(list(posterior = posterior, priorMLE = priorMLE))
}
