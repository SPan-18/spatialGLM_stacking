
sptGLM_exact <- function(y, X, S, time,
                         N.samp,
                         family = "poisson",
                         n_binom = NULL,
                         spCov = "matern",
                         mod_params,
                         Rfastparallel = FALSE,
                         verbose = TRUE){
  
  p <- dim(X)[2]
  n <- dim(X)[1]
  
  phi <- mod_params$phi
  phi_t <- mod_params$phi_t
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  nu_matern <- mod_params$nu_matern
  
  distS <- as.matrix(dist(S))
  distT <- as.matrix(dist(time))
  V_z <- matern(u = distS, phi = 1/phi, kappa = nu_matern) * exp(- phi_t * distT)
  # V_z <- diag(n)
  L_z <- chol(V_z)
  
  XtX <- crossprod(X)
  G <- makeG(X)
  XtG <- crossprod(X, G)
  if(verbose) cat("Big Cholesky..")
  HtH <- rbind(cbind(2*diag(n), X, G),
               cbind(t(X), XtX + diag(p), XtG),
               cbind(t(G), t(XtG), GtGplusI(X)))
  HtHchol <- chol(HtH)
  if(verbose) cat("done!\n")
  
  xi_samps <- array(dim = c(N.samp, n))
  beta_samps <- array(dim = c(N.samp, p))
  z_samps <- array(dim = c(N.samp, n*p))
  
  # progress bar from "progress" package
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = N.samp,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  for(i in 1:N.samp){
    
    # Updates the current state
    pb$tick()
    
    sigmasq_xi <- 1
    sigmasq_beta <- 1/rgamma(1, 0.5*nu_beta, 0.5*nu_beta)
    sigmasq_z <- 1/rgamma(1, 0.5*nu_z, 0.5*nu_z)
    
    w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                       alpha = y[x] + alpha_epsilon, 
                                       kappa = 1, psi = "psi3"))
    w_xi <- rnorm(n = n, mean = 0, sd = sqrt(sigmasq_xi))
    w_beta <- rnorm(n = p, mean = 0, sd = sqrt(sigmasq_beta))
    
    w_z <- rnorm(n = n*p, mean = 0, sd = sqrt(sigmasq_z))
    for(j in 1:p){
      ids <- (j-1)*n+1:n
      w_z[ids] <- as.numeric(crossprod(L_z, w_z[ids]))
    }
    
    Htw <- c(w_y + w_xi,
             as.numeric(crossprod(X, w_y)) + w_beta,
             as.numeric(crossprod(G, w_y)) + w_z)
    # # Htw <- crossprod(H, c(w_y, w_xi, w_beta, w_z))
    temp <- backsolve(HtHchol, Htw, transpose = TRUE)
    gamma <- backsolve(HtHchol, temp)
    
    # gamma <- projection(X = bigX, prechol = XtXplusIchol,
    #                     wy = w_y, wxi = w_xi, wbeta = w_beta, wz = w_z)
    
    xi_samps[i, ] <- gamma[1:n]
    beta_samps[i, ] <- gamma[(n+1):(n+p)]
    z_samps[i, ] <- gamma[(n+p+1):(n+p+n*p)]
  }
  return(list(beta = beta_samps,
              z = z_samps,
              xi = xi_samps))
}



