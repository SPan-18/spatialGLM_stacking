
# adaptive MCMC updates of spatiotemporal process parameters

sptParams_adaMCMC <- function(n, r, z, theta_starting, 
                              distS, distT,
                              n.batch, batch.length, n.report,
                              phi_s_Unifa, phi_s_Unifb, 
                              phi_t_Unifa, phi_t_Unifb, nu_z){
  
  target <- function(theta){
    phi_s <- logitinv_bounds(theta[1:r], phi_s_Unifa, phi_s_Unifb)
    phi_t <- logitinv_bounds(theta[(r+1):(2*r)], phi_t_Unifa, phi_t_Unifb)
    
    V_z <- lapply(1:r, function(x){ 1/(1 + phi_t[x]*distT) * 
        exp(- (phi_s[x]*distS) / sqrt(1 + phi_t[x]*distT)) })
    L_z <- lapply(V_z, 'chol')
    Q_z <- lapply(1:r, function(x){ 
      0.5*(nu_z + n)*log(1 + sum(backsolve(L_z[[x]], z[((x-1)*n+1):(x*n)], transpose = TRUE)^2)/nu_z) })
    
    logdetz <- lapply(L_z, function(x){ sum(log(diag(x))) })
    sumlogdetz <- do.call('sum', logdetz)
    sumlogdens <- do.call('sum', Q_z)
    logPost <- - sumlogdetz - sumlogdens
    logPost <- logPost + sum(log(phi_s - phi_s_Unifa)) + sum(log(phi_s_Unifb - phi_s))
    logPost <- logPost + sum(log(phi_t - phi_t_Unifa)) + sum(log(phi_t_Unifb - phi_t))
    
    return(logPost)
  }
  
  inits <- array(dim = 2*r)
  inits[1:r] <- logit_bounds(theta_starting[1:r], phi_s_Unifa, phi_s_Unifb)
  inits[(r+1):(2*r)] <- logit_bounds(theta_starting[(r+1):(2*r)], phi_t_Unifa, phi_t_Unifb)
  
  chain <- adaptMetropGibbs(ltd = target, starting = inits,
                            batch = n.batch, 
                            batch.length = batch.length,
                            verbose = FALSE)
  
  post_phi_s <- apply(chain$p.theta.samples[, 1:r], 1, 
                      function(x){logitinv_bounds(x, phi_s_Unifa, phi_s_Unifb)})
  post_phi_t <- apply(chain$p.theta.samples[, (r+1):(2*r)], 1, 
                      function(x){logitinv_bounds(x, phi_t_Unifa, phi_t_Unifb)})
  rownames(post_phi_s) <- NULL
  rownames(post_phi_t) <- NULL
  
  return(list(phi_s = post_phi_s,
              phi_t = post_phi_t))
}

# adaptive MCMC-within-Gibbs main function

sptvGLM_adaMetropGibbs <- function(y, X, X_tilde, S, time, 
                                   family, n_binom = NULL, 
                                   N.samp, starting, prior,
                                   n.batch = 3, batch.length = 10){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  r <- dim(X_tilde)[2]
  
  # priors for spatial process parameters
  phi_s_Unifa <- prior$phi_s_a
  phi_s_Unifb <- prior$phi_s_b
  phi_t_Unifa <- prior$phi_t_a
  phi_t_Unifb <- prior$phi_t_b
  
  # priors for variance parameters
  nu_xi <- prior$nu_xi
  nu_beta <- prior$nu_beta
  nu_z <- prior$nu_z
  
  # boundary adjustment parameter
  epsilon <- prior$alpha_epsilon
  
  nsp_params <- 2 * r
  
  post_phi_s <- array(dim = c(r, N.samp + 1))
  post_phi_t <- array(dim = c(r, N.samp + 1))
  post_xi <- array(dim = c(n, N.samp + 1))
  post_beta <- array(dim = c(p, N.samp + 1))
  post_z <- array(dim = c(n*r, N.samp + 1))
  
  post_phi_s[, 1] <- starting$phi_s
  post_phi_t[, 1] <- starting$phi_t
  
  post_xi[, 1] <- rep(0, n)
  post_beta[, 1] <- starting$beta
  post_z[, 1] <- rep(0, n*r)
  
  # pre-allocation of quantities used in projection
  diag1 <- apply(X_tilde, 1, function(x){ 1/(2 + sum(x^2)) })
  S_D_star <- crossprod(X, diag1 * X)
  diag(S_D_star) <- diag(S_D_star) + 1
  chol_p <- chol(S_D_star)
  chol_nr <- lapply(1:n, function(x){chol(tcrossprod(X_tilde[x, ]) + diag(r))})
  
  # Cholesky of first V_z
  distS <- dist(S)
  distT <- dist(time)
  distS <- dist2mat(distS, 128)
  distT <- dist2mat(distT, 128)
  V_z <- lapply(1:r, function(x){ 1/(1 + post_phi_t[x, 1]*distT) * 
      exp(- (post_phi_s[x, 1]*distS) / sqrt(1 + post_phi_t[x, 1]*distT)) })
  L_z <- lapply(V_z, 'chol')
  
  # n.batch <- 3
  # batch.length <- 10
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = N.samp,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = TRUE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  for(s in 1:N.samp){
    
    pb$tick()
    
    # sample posterior of eta and xi, beta, z from marginal priors
    if(family == "poisson"){
      v_eta <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                           alpha = y[x] + epsilon, 
                                           kappa = 1, psi = "psi3"))
    }else if(family == "binomial"){
      if(is.null(n_binom)){stop("Supply n_binom.")}
      v_eta <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                           alpha = y[x] + epsilon, 
                                           kappa = n_binom[x] + 2*epsilon, 
                                           psi = "psi2"))
    }else{ stop("Incorrect family.") }
    
    sigmasq_beta <- 1/rgamma(1, 0.5 * nu_beta, 0.5 * nu_beta)
    v_xi <- rnorm(n, mean = 0, sd = nu_xi)
    v_beta <- rnorm(p) * sqrt(sigmasq_beta)
    v_z <- sapply(1:r, function(x){ 
      sigmasq_z <- 1/rgamma(1, 0.5 * nu_z, 0.5 * nu_z)
      temp <- rnorm(n, mean = 0, sd = sqrt(sigmasq_z))
      temp <- crossprod(L_z[[x]], temp)
      return(temp)})
    v_z <- as.numeric(t(v_z))
    
    # evaluate projection
    gamma <- projection2(X = X, X_tilde = X_tilde,
                         chol_p = chol_p, chol_nr = chol_nr,
                         v_eta = v_eta, v_xi = v_xi,
                         v_beta = v_beta, v_z = v_z)
    temp_z <- matrix(gamma[(n+p+1):(n+p+n*r)], nrow = n, ncol = r, byrow = T)
    gamma[(n+p+1):(n+p+n*r)] <- as.numeric(c(temp_z))
    
    post_xi[, s+1] <- gamma[1:n]
    post_beta[, s+1] <- gamma[(n+1):(n+p)]
    post_z[, s+1] <- gamma[(n+p+1):(n+p+n*r)]
    
    # draw posterior of spatial process parameters
    
    theta_chain <- sptParams_adaMCMC(n = n, r = r, z = post_z[, s+1], 
                                     theta_starting = c(post_phi_s[, s], post_phi_t[, s]), 
                                     distS = distS, distT = distT, 
                                     n.batch = n.batch, 
                                     batch.length = batch.length,
                                     phi_s_Unifa = phi_s_Unifa, 
                                     phi_s_Unifb = phi_s_Unifb, 
                                     phi_t_Unifa = phi_t_Unifa, 
                                     phi_t_Unifb = phi_t_Unifb, 
                                     nu_z = nu_z)
    
    # id <- sample(n.batch*batch.length, 1)
    id <- n.batch*batch.length
    post_phi_s[, s+1] <- theta_chain$phi_s[, id]
    post_phi_t[, s+1] <- theta_chain$phi_t[, id]
  }
  
  return(list(beta = post_beta,
              z = post_z,
              xi = post_xi,
              phi_s = post_phi_s,
              phi_t = post_phi_t))
}
