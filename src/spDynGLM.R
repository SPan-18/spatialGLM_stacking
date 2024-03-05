
spDynGLM_stack <- function(y, Xt, S, time,
                           N.samp, MC.samp = 200,
                           family = "poisson",
                           n_binom = NULL,
                           spCov = "matern",
                           beta_prior = "nonstationary",
                           omega = NULL,
                           mod_params_list,
                           CV_fold = 10,
                           Rfastparallel = FALSE,
                           mc.cores = NULL,
                           verbose = TRUE,
                           solver = "ECOS"){
  
  distmat <- Rfast::Dist(S)
  
  n <- length(y)
  permut <- sample(1:n)
  y <- y[permut]
  time <- time[permut]
  if(!is.null(n_binom)){ n_binom <- n_binom[permut] }
  Xt <- Xt[permut, ]
  S <- S[permut, ]
  
  t_start <- Sys.time()
  
  if(is.null(mc.cores)){
    samps <- lapply(1:length(mod_params_list), function(x)
      spDynGLM_fixed(y = y, Xt = Xt, time = time,
                     distmat = distmat,
                     spCov = spCov,
                     N.samp = N.samp, MC.samp = MC.samp,
                     family = family,
                     n_binom = n_binom,
                     beta_prior = beta_prior,
                     omega = omega,
                     mod_params = mod_params_list[[x]],
                     CV_K = CV_fold, 
                     Rfastparallel = Rfastparallel))
  }else{
    
    ncores <- parallel::detectCores()
    if(mc.cores > ncores) warning("Number of requested cores exceeds limit (", ncores,"). 
                                  Setting mc.cores = ", ncores, ".")
    
    samps <- mclapply(1:length(mod_params_list), function(x)
      spDynGLM_fixed(y = y, Xt = Xt, time = time,
                     distmat = distmat,
                     spCov = spCov,
                     N.samp = N.samp, MC.samp = MC.samp,
                     family = family,
                     n_binom = n_binom,
                     beta_prior = beta_prior,
                     omega = omega,
                     mod_params = mod_params_list[[x]],
                     CV_K = CV_fold, 
                     Rfastparallel = Rfastparallel),
      mc.cores = min(mc.cores, ncores))
  }    # end multiple runs of models
  
  # for(i in 1:length(samps)){
  #   samps[[i]]$z <- samps[[i]]$z[order(permut), ]
  #   samps[[i]]$xi <- samps[[i]]$xi[order(permut), ]
  #   samps[[i]]$elpd <- samps[[i]]$elpd[order(permut)]
  # }
  # 
  # return(list(models = samps))
  
  elpd_mat <- do.call(cbind, lapply(samps, function(x) x$elpd))

  w_hat <- CVXR_stacking_weights(elpd_mat, solver = solver)
  w_hat <- as.numeric(w_hat)
  if(solver == "MOSEK"){
    w_hat <- sapply(w_hat, function(x) max(0, x))
    w_hat <- w_hat / sum(w_hat)
  }

  t_end <- Sys.time()

  runtime <- difftime(t_end, t_start)
  if(verbose) cat("\nRUNTIME:", round(runtime, 2), units(runtime), ".\n\n")

  if(verbose){
    stack_out <- as.matrix(do.call(rbind, lapply(mod_params_list, unlist)))
    stack_out <- cbind(stack_out, round(w_hat, 2))
    colnames(stack_out) = c("phi", "smooth", "phi_t", "rho", "epsilon", "nu.xi", "nu.beta", "nu.z", "weight")
    rownames(stack_out) = paste("Model", 1:nrow(stack_out))
    cat("STACKING WEIGHTS:\n")
    print(knitr::kable(stack_out))
  }

  for(i in 1:length(samps)){
    samps[[i]]$z <- samps[[i]]$z[order(permut), ]
    samps[[i]]$xi <- samps[[i]]$xi[order(permut), ]
    samps[[i]]$elpd <- samps[[i]]$elpd[order(permut)]
  }

  return(list(models = samps, weights = w_hat))
  
}    # end spDynGLM function



spDynGLM_fixed <- function(y, Xt, time,
                           distmat,
                           N.samp, MC.samp,
                           spCov = "matern",
                           n_binom = NULL,
                           family, 
                           beta_prior = "nonstationary",
                           omega = NULL,
                           mod_params,
                           CV_K = 10,
                           Rfastparallel = FALSE){
  
  n <- dim(Xt)[1]
  p <- dim(Xt)[2]
  Nt <- max(time)
  
  phi <- mod_params$phi
  phi_t <- mod_params$phi_t
  rho <- mod_params$rho
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  nu_matern <- mod_params$nu_matern
  if(beta_prior == "nonstationary"){
    omega <- omega
  }
  
  distT <- as.matrix(dist(time + rnorm(n, 0, 0.01)))
  
  if(spCov == "exponential"){
    V_z <- exp(- phi * distmat) * exp(- phi_t * Rfast::Dist(time))
  }else if(spCov == "matern"){
    V_z <- matern(u = distmat, phi = 1/phi, kappa = nu_matern) * exp(- phi_t * distT)
    # V_z <- matern(u = distmat, phi = 1/phi, kappa = nu_matern) * phi_t^as.matrix(dist(time))
  }else{
    V_z <- geoR::cov.spatial(distmat, cov.model = spCov, 
                             cov.pars = c(1, 1/phi), kappa = nu_matern) * exp(- phi_t * Rfast::Dist(time))
  }
  
  L_z <- Rfast::cholesky(V_z, parallel = Rfastparallel)
  X <- mkdynX(Xt, time)
  XtXplusI <- crossprod(X) / 3 + diag(p * Nt)
  XtXplusIchol <- chol(XtXplusI)
  
  gamma <- sapply(1:N.samp, function(x) 
    sampler_dynGCM(Nt = Nt, y = y, X = X, 
                   family = family, n_binom = n_binom, 
                   XtXplusIchol = XtXplusIchol, 
                   L_z = L_z, 
                   beta_prior = beta_prior,
                   rho = rho, omega = omega,
                   nu_xi = nu_xi, nu_beta = nu_beta, nu_z = nu_z, 
                   alpha_epsilon = alpha_epsilon))
  
  # return(list(beta = gamma[(n+1):(n+(p*Nt)), ],
  #             z = gamma[(n+(p*Nt)+1):(2*n+(p*Nt)), ],
  #             xi = gamma[1:n, ]))
  
  elpd <- CV_dynsampler(Nt = Nt, y = y, X = X,
                        N.samp = MC.samp,
                        V_z_full = V_z, L_z_full = L_z,
                        family = family,
                        n_binom = n_binom,
                        beta_prior = beta_prior,
                        omega = omega,
                        mod_params = mod_params,
                        CV_K = CV_K,
                        Rfastparallel = Rfastparallel)

  return(list(beta = gamma[(n+1):(n+(p*Nt)), ],
              z = gamma[(n+(p*Nt)+1):(2*n+(p*Nt)), ],
              xi = gamma[1:n, ],
              elpd = elpd))
  
}

sampler_dynGCM <- function(Nt, y, X, family,
                           n_binom = NULL,
                           XtXplusIchol, L_z,
                           beta_prior = "nonstationary",
                           rho, omega,
                           nu_xi, nu_beta, nu_z, alpha_epsilon){
  
  n <- length(y)
  p <- dim(X)[2] / Nt
  
  w_xi <- rnorm(n, mean = 0, sd = nu_xi)
  # sigmasq_beta <- 1/rgamma(1, 0.5 * nu_beta, 0.5 * nu_beta)
  sigmasq_z <- 1/rgamma(1, 0.5 * nu_z, 0.5 * nu_z)
  
  if(beta_prior == "stationary"){
    sigmasq_beta <- 1/rgamma(1, 0.5 * nu_beta, 0.5 * nu_beta)
    w_beta <- rnorm(p*Nt) * sqrt(sigmasq_beta)
    L_beta <- mkbetacov(size = Nt, rho = rho)
    L_beta <- chol(L_beta)
    for(i in 1:p){
      ids <- ((i-1)*Nt+1):(i*Nt)
      w_beta[ids] <- as.numeric(crossprod(L_beta, w_beta[ids]))
    }
  }else if(beta_prior == "nonstationary"){
    
    vt <- array(dim = Nt)
    nt <- array(dim = Nt)
    for(t in 1:Nt){
      if(t == 1){
        nt[t] <- nu_beta
        vt[t] <- 1 * omega / rbeta(1, omega*nt[t], (1 - omega)*nt[t])
      }else{
        nt[t] <- nt[t-1] + 1
        vt[t] <- vt[t-1] * omega / rbeta(1, omega*nt[t], (1 - omega)*nt[t])
      }
    }
    sigmasq_beta <- rep(vt, times = p)
    w_beta <- rnorm(p*Nt) * sqrt(sigmasq_beta)
    L_beta <- mkbetacov(size = Nt, rho = rho)
    L_beta <- chol(L_beta)
    for(i in 1:p){
      ids <- ((i-1)*Nt+1):(i*Nt)
      w_beta[ids] <- as.numeric(crossprod(L_beta, w_beta[ids]))
    }
  }
  
  if(family == "poisson"){
    # browser()
    w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                       alpha = y[x] + alpha_epsilon, 
                                       kappa = 1, psi = "psi3"))
  }else if(family == "binomial"){
    if(is.null(n_binom)){stop("Supply n_binom.")}
    w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                       alpha = y[x] + alpha_epsilon, 
                                       kappa = n_binom[x] + 2*alpha_epsilon, 
                                       psi = "psi2"))
  }else{ stop("Incorrect family.") }
  
  w_z <- rnorm(n) * sqrt(sigmasq_z)
  w_z <- as.numeric(crossprod(L_z, w_z))
  
  gamma <- projection(X = X, prechol = XtXplusIchol,
                      wy = w_y, wxi = w_xi, wbeta = w_beta, wz = w_z)
  
  return(gamma)
}


CV_dynsampler <- function(Nt, y, X,
                          N.samp,
                          V_z_full, L_z_full = NULL,
                          family,
                          n_binom = NULL,
                          beta_prior = "nonstationary",
                          omega,
                          mod_params,
                          CV_K = 10,
                          Rfastparallel = FALSE){
  
  n <- length(y)
  partition_list <- id_partition(n, CV_K, random = FALSE)
  
  CV_samps <- lapply(1:length(partition_list), function(x)
    elpd_dynGCM(Nt = Nt, 
                y_train = y[-partition_list[[x]]],
                X_train = X[-partition_list[[x]], ],
                y_pred = y[partition_list[[x]]],
                X_pred = X[partition_list[[x]], ],
                J_tilde = V_z_full[- partition_list[[x]],
                                   partition_list[[x]]],
                V_tilde = V_z_full[partition_list[[x]],
                                   partition_list[[x]]],
                V_z_train = V_z_full[-partition_list[[x]],
                                     -partition_list[[x]]],
                L_z_train = cholesky_CV(L_z_full, partition_list[[x]]),
                N.samp = N.samp,
                mod_params = mod_params,
                family = family,
                n_binom_train = n_binom[-partition_list[[x]]],
                n_binom_pred = n_binom[partition_list[[x]]],
                beta_prior = beta_prior,
                omega = omega,
                Rfastparallel = Rfastparallel))
  
  elpd <- array(dim = n)
  for(k in 1:CV_K){
    elpd[partition_list[[k]]] = CV_samps[[k]]
  }
  return(elpd)
}


elpd_dynGCM <- function(Nt, 
                        y_train, X_train, y_pred, X_pred, N.samp,
                        J_tilde, V_tilde, V_z_train, 
                        L_z_train = NULL,
                        family,
                        n_binom_train = NULL,
                        n_binom_pred = NULL,
                        beta_prior = "nonstationary",
                        omega = NULL,
                        mod_params,
                        Rfastparallel = FALSE){
  
  n <- length(y_train)
  p <- dim(X_train)[2] / Nt
  
  rho <- mod_params$rho
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  nu_matern <- mod_params$nu_matern
  if(beta_prior == "nonstationary"){
    omega <- omega
  }
  
  # attack here
  if(is.null(L_z_train)) L_z_train <- Rfast::cholesky(V_z_train, parallel = Rfastparallel)
  XtXplusI <- crossprod(X_train) / 3 + diag(p*Nt)
  XtXplusIchol <- chol(XtXplusI)
  
  gamma_train <- sapply(1:N.samp, function(x)  
    sampler_dynGCM(Nt = Nt, y = y_train, X = X_train, 
                   family = family,
                   n_binom = n_binom_train,
                   XtXplusIchol = XtXplusIchol, 
                   L_z = L_z_train,
                   beta_prior = beta_prior,
                   rho = rho, omega = omega,
                   nu_xi = nu_xi, 
                   nu_beta = nu_beta, 
                   nu_z = nu_z, 
                   alpha_epsilon = alpha_epsilon))
  
  # samp_xi <- gamma_train[1:n, ]
  samp_beta <- gamma_train[(n+1):(n+p*Nt), ]
  samp_z <- gamma_train[(n+p*Nt+1):(2*n+p*Nt), ]
  
  z_tilde <- predict_z(z_post = samp_z, 
                       J = J_tilde, cholV = L_z_train,
                       V_tilde = V_tilde, nu_z = nu_z)
  mu <- exp(X_pred %*% samp_beta + z_tilde)
  if(family == "poisson"){
    elpd_mat <- dpois(y_pred, mu, log = FALSE)
  }else if(family == "binomial"){
    prob <- ilogit(mu)
    elpd_mat <- dbinom(y_pred, n_binom_pred, prob, log = FALSE)
  }
  
  elpd <- apply(elpd_mat, 1, function(x) mean(x, na.rm = FALSE))
  
  return(log(elpd))
}

spDynGLM_exact <- function(y, X, S, time,
                           N.samp,
                           family = "poisson",
                           n_binom = NULL,
                           spCov = "matern",
                           beta_prior = "nonstationary",
                           omega = NULL,
                           mod_params,
                           Rfastparallel = FALSE,
                           verbose = TRUE){

  p <- dim(X)[2]
  n <- dim(X)[1]
  Nt <- max(time)
  
  phi <- mod_params$phi
  phi_t <- mod_params$phi_t
  rho <- mod_params$rho
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  nu_matern <- mod_params$nu_matern
  if(beta_prior == "nonstationary"){
    omega <- omega
  }

  distS <- Rfast::Dist(S)
  distT <- Rfast::Dist(time)
  V_z <- matern(u = distS, phi = 1/phi, kappa = nu_matern) * exp(- phi_t * distT)
  # V_z <- diag(n)
  L_z <- Rfast::cholesky(V_z)
  # L_beta <- rho^as.matrix(dist(1:Nt))
  L_beta <- mkbetacov(size = Nt, rho = rho)
  L_beta <- chol(L_beta)
  
  bigX <- mkdynX(X, time)
  # bigX <- mkdynX2(X, time)
  XtX <- crossprod(bigX)
  XtXplusI <- XtX / 3 + diag(p * Nt)
  XtXplusIchol <- chol(XtXplusI)
  # H <- rbind(cbind(diag(n), bigX, diag(n)),
  #            cbind(diag(n), matrix(0,n,p*Nt), matrix(0,n,n)),
  #            cbind(matrix(0,p*Nt,n), diag(p*Nt), matrix(0,p*Nt,n)),
  #            cbind(matrix(0,n,n), matrix(0,n,p*Nt), diag(n)))
  # HtH <- crossprod(H)
  # HtH <- rbind(cbind(2*diag(n), bigX, diag(n)),
  #              cbind(t(bigX), XtX + diag(p*Nt), t(bigX)),
  #              cbind(diag(n), bigX, 2*diag(n)))
  # HtHchol <- chol(HtH)
  
  beta_samps <- array(dim = c(N.samp, p*Nt))
  z_samps <- array(dim = c(N.samp, n))
  xi_samps <- array(dim = c(N.samp, n))
  
  for(i in 1:N.samp){
    sigmasq_xi <- 1
    # sigmasq_beta <- 1/rgamma(1, 0.5*nu_beta, 0.5*nu_beta)
    sigmasq_z <- 1/rgamma(1, 0.5*nu_z, 0.5*nu_z)
    
    w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                       alpha = y[x] + alpha_epsilon, 
                                       kappa = 1, psi = "psi3"))
    w_xi <- rnorm(n = n, mean = 0, sd = sqrt(sigmasq_xi))
    w_beta <- rnorm(n = p*Nt, mean = 0, sd = 1)
    
    if(beta_prior == "nonstationary"){
      vt <- array(dim = Nt)
      nt <- array(dim = Nt)
      for(t in 1:Nt){
        if(t == 1){
          nt[t] <- nu_beta
          vt[t] <- 1 * omega / rbeta(1, omega*nt[t], (1 - omega)*nt[t])
        }else{
          nt[t] <- nt[t-1] + 1
          vt[t] <- vt[t-1] * omega / rbeta(1, omega*nt[t], (1 - omega)*nt[t])
        }
      }
      sigmasq_beta <- rep(vt, times = p)
    }
    w_beta <- w_beta * sqrt(sigmasq_beta)
    # for(j in 1:p){
    #   id_t <- ((j-1)*Nt+1):(j*Nt)
    #   w_beta[id_t] <- as.numeric(crossprod(L_beta, w_beta[id_t]))
    #   w_beta[id_t] <- w_beta[id_t] * sqrt(sigmasq_beta)
    # }
    w_z <- rnorm(n = n, mean = 0, sd = sqrt(sigmasq_z))
    w_z <- as.numeric(crossprod(L_z, w_z))
    
    # Htw <- c(w_y + w_xi,
    #          as.numeric(crossprod(bigX, w_y)) + w_beta,
    #          w_y + w_z)
    # # Htw <- crossprod(H, c(w_y, w_xi, w_beta, w_z))
    # temp <- backsolve(HtHchol, Htw, transpose = TRUE)
    # gamma <- backsolve(HtHchol, temp)
    
    gamma <- projection(X = bigX, prechol = XtXplusIchol,
                        wy = w_y, wxi = w_xi, wbeta = w_beta, wz = w_z)
    
    xi_samps[i, ] <- gamma[1:n]
    beta_samps[i, ] <- gamma[(n+1):(n+p*Nt)]
    z_samps[i, ] <- gamma[(n+p*Nt+1):(2*n+p*Nt)]
  }
  return(list(beta = beta_samps,
              z = z_samps,
              xi = xi_samps))
}

