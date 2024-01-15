logit_bounds <- function(x, a, b){
  return(log(x - a) - log(b - x))
}

logitinv_bounds <- function(x, a, b){
  return(a + (b - a) * (1 / (1 + exp(- x))))
  # return(b - (b - a) / (1 + exp(x)))
}

spGCM_MCMC <- function(y, X, S, N.samp, family, 
                       spCov,
                       n_binom = NULL,
                       marginalPrior = TRUE,
                       starting,
                       prior,
                       tuning,
                       Rfastparallel = FALSE){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # prior for spatial process parameters
  phiUnifa <- prior$phi[1]
  phiUnifb <- prior$phi[2]
  nuUnifa <- prior$nu[1]
  nuUnifb <- prior$nu[2]
  
  # prior for variance parameters
  nu_xi <- prior$nu_xi
  nu_beta <- prior$nu_beta
  nu_z <- prior$nu_z
  
  # boundary adjustment parameter
  alpha_epsilon <- prior$alpha_epsilon
  
  if(marginalPrior){
    if(spCov == "exponential"){
      nsp_params <- 1
    }else if(spCov == "matern"){
      nsp_params <- 2
    }else{ stop("Correlation kernel not supported.")}
  }
  
  sp_params <- array(dim = c(nsp_params, N.samp + 1))
  post_phi <- array(dim = N.samp + 1)
  post_nu <- array(dim = N.samp + 1)
  post_xi <- array(dim = c(n, N.samp + 1))
  post_beta <- array(dim = c(p, N.samp + 1))
  post_z <- array(dim = c(n, N.samp + 1))
  
  post_phi[1] <- starting$phi
  if(spCov == "matern"){
    post_nu[1] <- starting$nu
  }
  post_xi[, 1] <- rep(0, n)
  post_beta[, 1] <- starting$beta
  post_z[, 1] <- rep(0, n)
  
  sp_params[1, 1] <- logit_bounds(post_phi[1], phiUnifa, phiUnifb)
  if(spCov == "matern"){
    sp_params[2, 1] <- logit_bounds(post_nu[1], nuUnifa, nuUnifb)
  }
  
  # Cholesky pre-processing
  distmat <- Rfast::Dist(S)
  if(spCov == "exponential"){
    V_z <- matern(u = distmat, phi = 1/post_phi[1], kappa = 0.5)
  }else if(spCov == "matern"){
    V_z <- matern(u = distmat, phi = 1/post_phi[1], kappa = post_nu[1])
  }else{
    stop("Correlation kernel not supported.")
  }
  L_z <- Rfast::cholesky(V_z, parallel = Rfastparallel)
  XtXplusI <- crossprod(X) / 3 + diag(p)
  XtXplusIchol <- chol(XtXplusI)
  
  logPostcurrent <- -Inf
  nAccept <- 0
  
  cat("Tuning: ", tuning, "\n")
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = N.samp,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  for(s in 1:N.samp){
    
    pb$tick()
    
    # draw posterior of (gamma,q) by Exact Posterior Regression
    
    w_xi <- rnorm(n, mean = 0, sd = nu_xi)
    sigmasq_beta <- 1/rgamma(1, 0.5 * nu_beta, 0.5 * nu_beta)
    sigmasq_z <- 1/rgamma(1, 0.5 * nu_z, 0.5 * nu_z)
    
    if(family == "poisson"){
      w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                         alpha = y[x] + alpha_epsilon, 
                                         kappa = 1, psi = "psi3"))
    }else if(family == "binomial"){
      if(is.null(n_binom)){stop("Supply n_binom.")}
      w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                         alpha = y[x] + alpha_epsilon, 
                                         kappa = n_binom[x] + 2*alpha_epsilon, 
                                         psi = "psi2"))
    }else{ stop("Family currently not supported.") }
    
    w_beta <- rnorm(p) * sqrt(sigmasq_beta)
    w_z <- rnorm(n) * sqrt(sigmasq_z)
    w_z <- as.numeric(crossprod(L_z, w_z))
    
    gamma <- projection(X = X, prechol = XtXplusIchol, 
                        wy = w_y, wxi = w_xi, 
                        wbeta = w_beta, wz = w_z)
    
    post_xi[, s+1] <- gamma[1:n, ]
    post_beta[, s+1] <- gamma[(n+1):(n+p), ]
    post_z[, s+1] <- gamma[(n+p+1):(2*n+p), ]
    
    
    # draw posterior of spatial process parameters
    
    sp_params_cand <- array(dim = nsp_params)
    for(i in 1:nsp_params){
      sp_params_cand[i] <- rnorm(1, sp_params[i, s], tuning[i])  # needs tuning 
    }
    phi <- logitinv_bounds(sp_params_cand[1], phiUnifa, phiUnifb)
    if(spCov == "matern"){
      nu <- logitinv_bounds(sp_params_cand[2], nuUnifa, nuUnifb)
    }
    
    if(spCov == "exponential"){
      V_z_cand <- matern(u = distmat, phi = 1/phi, kappa = 0.5)
    }else if(spCov == "matern"){
      V_z_cand <- matern(u = distmat, phi = 1/phi, kappa = nu)
    }
    L_z_cand <- Rfast::cholesky(V_z_cand, parallel = Rfastparallel)
    
    # evaluate log posterior of candidate phi and nu
    
    logPostcand <- 0.0
    logdetz <- 2 * sum(log(diag(L_z_cand)))
    Qz <- sum(backsolve(L_z_cand, post_z[, s+1], transpose = TRUE)^2)
    logPostcand <- logPostcand - 0.5 * logdetz - 
      0.5 * (nu_z + n) * log(1 + Qz / nu_z)
    logPostcand <- logPostcand + log(phi - phiUnifa) + log(phiUnifb - phi)
    if(spCov == "matern"){
      logPostcand <- logPostcand + log(nu - nuUnifa) + log(nuUnifb - nu)
    }
    
    # calculate Metropolis-Hastings ratio
    
    logMHratio <- logPostcand - logPostcurrent
    
    if(runif(1, 0.0, 1.0) <= exp(logMHratio)){
      sp_params[, s+1] <- sp_params_cand
      L_z <- L_z_cand
      logPostcurrent <- logPostcand
      nAccept <- nAccept + 1
    }else{
      sp_params[, s+1] <- sp_params[, s]
    }
    
    post_phi[s+1] <- logitinv_bounds(sp_params[1, s+1], phiUnifa, phiUnifb)
    if(spCov == "matern"){
      post_nu[s+1] <- logitinv_bounds(sp_params[2, s+1], nuUnifa, nuUnifb)
    }
  }
  
  cat("Acceptance = ", nAccept / N.samp * 100, "%\n")
  
  return(list(beta = post_beta,
              z = post_z,
              xi = post_xi,
              phi = post_phi,
              nu = post_nu))
}



spGCM_MCMC2 <- function(y, X, S, N.samp, family, 
                        spCov,
                        n_binom = NULL,
                        marginalPrior = TRUE,
                        starting,
                        prior,
                        tuning,
                        Rfastparallel = FALSE){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # prior for spatial process parameters
  phiUnifa <- prior$phi[1]
  phiUnifb <- prior$phi[2]
  nuUnifa <- prior$nu[1]
  nuUnifb <- prior$nu[2]
  
  # prior for variance parameters
  nu_xi <- prior$nu_xi
  nu_beta <- prior$nu_beta
  nu_z <- prior$nu_z
  
  # boundary adjustment parameter
  alpha_epsilon <- prior$alpha_epsilon
  
  if(marginalPrior){
    if(spCov == "exponential"){
      nsp_params <- 1
    }else if(spCov == "matern"){
      nsp_params <- 2
    }else{ stop("Correlation kernel not supported.")}
  }
  
  sp_params <- array(dim = c(nsp_params, N.samp + 1))
  sp_params_cand <- array(dim = nsp_params)
  post_phi <- array(dim = N.samp + 1)
  post_nu <- array(dim = N.samp + 1)
  post_xi <- array(dim = c(n, N.samp + 1))
  post_beta <- array(dim = c(p, N.samp + 1))
  post_z <- array(dim = c(n, N.samp + 1))
  
  post_phi[1] <- starting$phi
  if(spCov == "matern"){
    post_nu[1] <- starting$nu
  }
  post_xi[, 1] <- rep(0, n)
  post_beta[, 1] <- starting$beta
  post_z[, 1] <- rep(0, n)
  
  sp_params[1, 1] <- logit_bounds(post_phi[1], phiUnifa, phiUnifb)
  if(spCov == "matern"){
    sp_params[2, 1] <- logit_bounds(post_nu[1], nuUnifa, nuUnifb)
  }
  
  # Cholesky pre-processing
  distmat <- Rfast::Dist(S)
  if(spCov == "exponential"){
    V_z <- matern(u = distmat, phi = 1/post_phi[1], kappa = 0.5)
  }else if(spCov == "matern"){
    V_z <- matern(u = distmat, phi = 1/post_phi[1], kappa = post_nu[1])
  }else{
    stop("Correlation kernel not supported.")
  }
  L_z <- Rfast::cholesky(V_z, parallel = Rfastparallel)
  XtXplusI <- crossprod(X) / 3 + diag(p)
  XtXplusIchol <- chol(XtXplusI)
  
  logPostcurrent <- -Inf
  nAccept <- 0
  
  for(s in 1:N.samp){
    
    cat("Running iteration", s, "\n")
    
    # draw posterior of (gamma,q) by Exact Posterior Regression
    
    w_xi <- rnorm(n, mean = 0, sd = nu_xi)
    sigmasq_beta <- 1/rgamma(1, 0.5 * nu_beta, 0.5 * nu_beta)
    sigmasq_z <- 1/rgamma(1, 0.5 * nu_z, 0.5 * nu_z)
    
    if(family == "poisson"){
      w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                         alpha = y[x] + alpha_epsilon, 
                                         kappa = 1, psi = "psi3"))
    }else if(family == "binomial"){
      if(is.null(n_binom)){stop("Supply n_binom.")}
      w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                         alpha = y[x] + alpha_epsilon, 
                                         kappa = n_binom[x] + 2*alpha_epsilon, 
                                         psi = "psi2"))
    }else{ stop("Family currently not supported.") }
    
    w_beta <- rnorm(p) * sqrt(sigmasq_beta)
    w_z <- rnorm(n) * sqrt(sigmasq_z)
    w_z <- as.numeric(crossprod(L_z, w_z))
    
    gamma <- projection(X = X, prechol = XtXplusIchol, 
                        wy = w_y, wxi = w_xi, 
                        wbeta = w_beta, wz = w_z)
    
    post_xi[, s+1] <- gamma[1:n, ]
    post_beta[, s+1] <- gamma[(n+1):(n+p), ]
    post_z[, s+1] <- gamma[(n+p+1):(2*n+p), ]
    
    
    # draw posterior of spatial process parameters
    
    accept <- 0
    
    while(accept == 0){
      for(i in 1:nsp_params){
        sp_params_cand[i] <- rnorm(1, sp_params[i, s], tuning[i])  # needs tuning
      }
      phi <- logitinv_bounds(sp_params_cand[1], phiUnifa, phiUnifb)
      if(spCov == "matern"){
        nu <- logitinv_bounds(sp_params_cand[2], nuUnifa, nuUnifb)
      }

      if(spCov == "exponential"){
        V_z_cand <- matern(u = distmat, phi = 1/phi, kappa = 0.5)
      }else if(spCov == "matern"){
        V_z_cand <- matern(u = distmat, phi = 1/phi, kappa = nu)
      }
      L_z_cand <- Rfast::cholesky(V_z_cand, parallel = Rfastparallel)

      # evaluate log posterior of candidate phi and nu

      logPostcand <- 0.0
      logdetz <- 2 * sum(log(diag(L_z_cand)))
      Qz <- sum(backsolve(L_z_cand, post_z[, s+1], transpose = TRUE)^2)
      logPostcand <- logPostcand - 0.5 * logdetz -
        0.5 * (nu_z + n) * log(1 + Qz / nu_z)
      logPostcand <- logPostcand + log(phi - phiUnifa) + log(phiUnifb - phi)
      if(spCov == "matern"){
        logPostcand <- logPostcand + log(nu - nuUnifa) + log(nuUnifb - nu)
      }
      
      logMHratio <- logPostcand - logPostcurrent
      cat(phi, logPostcand, logPostcurrent, logMHratio, round(exp(logMHratio), 3), logdetz, log(1+Qz/nu_z))
      if(runif(1, 0.0, 1.0) <= exp(logMHratio)){
        cat("\t A \n")
        sp_params[, s+1] <- sp_params_cand
        L_z <- L_z_cand
        logPostcurrent <- logPostcand
        accept <- 1
      }else{
        cat("\t R \n")
      }
    }
    
    post_phi[s+1] <- logitinv_bounds(sp_params[1, s+1], phiUnifa, phiUnifb)
    if(spCov == "matern"){
      post_nu[s+1] <- logitinv_bounds(sp_params[2, s+1], nuUnifa, nuUnifb)
    }
  }
  
  return(list(beta = post_beta,
              z = post_z,
              xi = post_xi,
              phi = post_phi,
              nu = post_nu))
}
