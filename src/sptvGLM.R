
sptvGLM_exact <- function(y, X, X_tilde, S, time,
                          N.samp,
                          family = "poisson",
                          n_binom = NULL,
                          mod_params,
                          verbose = TRUE){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  r <- dim(X_tilde)[2]
  
  phi_s <- mod_params$phi_s
  phi_t <- mod_params$phi_t
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  
  t_start <- Sys.time()
  # evaluate correlation matrix of spatial-temporal process
  distS <- as.matrix(dist(S))
  distT <- as.matrix(dist(time))
  V_z <- 1/(1 + phi_t*distT) * exp(- (phi_s*distS) / sqrt(1 + phi_t*distT))
  L_z <- chol(V_z)
  
  # pre-allocation of useful Cholesky
  diag1 <- apply(X_tilde, 1, function(x){ 1/(2 + sum(x^2)) })
  S_D_star <- crossprod(X, diag1 * X)
  diag(S_D_star) <- diag(S_D_star) + 1
  chol_p <- chol(S_D_star)
  chol_nr <- lapply(1:n, function(x){chol(tcrossprod(X_tilde[x, ]) + diag(r))})
  
  # sample from the posterior
  gamma <- sapply(1:N.samp, function(x) 
    sampler_sptv(y = y, X = X, X_tilde = X_tilde,
                 family = family, n_binom = n_binom, 
                 chol_p = chol_p, chol_nr = chol_nr, L_z = L_z, 
                 nu_xi = nu_xi, nu_beta = nu_beta, nu_z = nu_z, 
                 alpha_epsilon = alpha_epsilon))
  t_end <- Sys.time()
  t_elapsed <- t_end - t_start
  if(verbose){
    cat("Time elapsed =", t_elapsed, units(t_elapsed), "\n")
  }
  
  return(list(xi = gamma[1:n, ],
              beta = gamma[(n+1):(n+p), ],
              z = gamma[(n+p+1):(n+p+n*r), ]))
}

sampler_sptv <- function(y, X, X_tilde, family, n_binom = NULL,
                         chol_p, chol_nr, L_z,
                         nu_xi, nu_beta, nu_z, alpha_epsilon){
  
  n <- length(y)
  p <- dim(X)[2]
  r <- dim(X_tilde)[2]
  
  # sample posterior of eta and xi, beta, z from marginal priors
  if(family == "poisson"){
    v_eta <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                         alpha = y[x] + alpha_epsilon, 
                                         kappa = 1, psi = "psi3"))
  }else if(family == "binomial"){
    if(is.null(n_binom)){stop("Supply n_binom.")}
    v_eta <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                         alpha = y[x] + alpha_epsilon, 
                                         kappa = n_binom[x] + 2*alpha_epsilon, 
                                         psi = "psi2"))
  }else{ stop("Incorrect family.") }
  
  sigmasq_beta <- 1/rgamma(1, 0.5 * nu_beta, 0.5 * nu_beta)
  # sigmasq_z <- 1/rgamma(1, 0.5 * nu_z, 0.5 * nu_z)
  v_xi <- rnorm(n, mean = 0, sd = nu_xi)
  v_beta <- rnorm(p) * sqrt(sigmasq_beta)
  v_z <- sapply(1:r, function(x){ 
    sigmasq_z <- 1/rgamma(1, 0.5 * nu_z, 0.5 * nu_z)
    temp <- rnorm(n, mean = 0, sd = sqrt(sigmasq_z))
    temp <- crossprod(L_z, temp)
    return(temp)})
  v_z <- as.numeric(t(v_z))
  
  # evaluate projection
  gamma <- projection2(X = X, X_tilde = X_tilde,
                      chol_p = chol_p, chol_nr = chol_nr,
                      v_eta = v_eta, v_xi = v_xi,
                      v_beta = v_beta, v_z = v_z)
  temp_z <- matrix(gamma[(n+p+1):(n+p+n*r)], nrow = n, ncol = r, byrow = T)
  gamma[(n+p+1):(n+p+n*r)] <- as.numeric(c(temp_z))
  
  return(gamma)
}


