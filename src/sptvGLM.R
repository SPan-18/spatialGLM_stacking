
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
  epsilon <- mod_params$epsilon
  
  t_start <- Sys.time()
  # evaluate correlation matrix of spatial-temporal process
  # distS <- as.matrix(dist(S))
  # distT <- as.matrix(dist(time))
  distS <- dist(S)
  distT <- dist(time)
  distS <- dist2mat(distS, 128)
  distT <- dist2mat(distT, 128)
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
                 epsilon = epsilon))
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
                         nu_xi, nu_beta, nu_z, epsilon){
  
  n <- length(y)
  p <- dim(X)[2]
  r <- dim(X_tilde)[2]
  
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

sptvGLM_stack <- function(y, X, X_tilde, S, time,
                          N.samp, MC.samp = 200,
                          family = "poisson", n_binom = NULL,
                          mod_params_list,
                          CV_fold = 10, mc.cores = NULL, 
                          solver = "ECOS", verbose = TRUE){
  
  n <- length(y)
  p <- dim(X)[2]
  r <- dim(X_tilde)[2]
  permut <- sample(1:n)
  y <- y[permut]
  time <- time[permut]
  if(!is.null(n_binom)){ n_binom <- n_binom[permut] }
  X <- X[permut, ]
  X_tilde <- X_tilde[permut, ]
  S <- S[permut, ]
  time <- time[permut]
  
  distS <- dist(S)
  distT <- dist(time)
  distS <- dist2mat(distS, 128)
  distT <- dist2mat(distT, 128)
  
  t_start <- Sys.time()
  
  if(is.null(mc.cores)){
    samps <- lapply(1:length(mod_params_list), function(x)
      sptvGLM_fixed(y = y, X = X, X_tilde = X_tilde,
                    distmatS = distS, distmatT = distT,
                    N.samp = N.samp, MC.samp = MC.samp,
                    family = family, n_binom = n_binom,
                    mod_params = mod_params_list[[x]],
                    CV_K = CV_fold))
  }else{
    
    ncores <- parallel::detectCores()
    if(mc.cores > ncores) warning("Number of requested cores exceeds limit (", ncores,"). 
                                  Setting mc.cores = ", ncores, ".")
    
    samps <- mclapply(1:length(mod_params_list), function(x)
      sptvGLM_fixed(y = y, X = X, X_tilde = X_tilde,
                    distmatS = distS, distmatT = distT,
                    N.samp = N.samp, MC.samp = MC.samp,
                    family = family, n_binom = n_binom,
                    mod_params = mod_params_list[[x]],
                    CV_K = CV_fold),
      mc.cores = min(mc.cores, ncores))
  }    # end multiple runs of models
  
  elpd_mat <- do.call(cbind, lapply(samps, function(x) x$elpd))
  
  # for(i in 1:length(samps)){
  #   samps[[i]]$z <- samps[[i]]$z[order(permut), ]
  #   samps[[i]]$xi <- samps[[i]]$xi[order(permut), ]
  #   samps[[i]]$elpd <- samps[[i]]$elpd[order(permut)]
  # }
  
  # return(elpd_mat)
  
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
    stack_out_names <- colnames(stack_out)
    stack_out <- cbind(stack_out, round(w_hat, 3))
    colnames(stack_out) = c(stack_out_names, "weight")
    rownames(stack_out) = paste("Model", 1:nrow(stack_out))
    cat("STACKING WEIGHTS:\n")
    print(knitr::kable(stack_out))
  }
  
  for(i in 1:length(samps)){
    for(j in 1:r){
      ids <- ((j-1)*n+1):(j*n)
      temp <- samps[[i]]$z[ids, ]
      temp <- temp[order(permut), ]
      samps[[i]]$z[ids, ] <- temp
    }
    samps[[i]]$xi <- samps[[i]]$xi[order(permut), ]
    samps[[i]]$elpd <- samps[[i]]$elpd[order(permut)]
  }

  return(list(models = samps, weights = w_hat))
}

sptvGLM_fixed <- function(y, X, X_tilde,
                          distmatS, distmatT,
                          N.samp, MC.samp,
                          family, n_binom,
                          mod_params, CV_K = 10){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  r <- dim(X_tilde)[2]
  
  phi_s <- mod_params$phi_s
  phi_t <- mod_params$phi_t
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  epsilon <- mod_params$epsilon
  
  V_z <- 1/(1 + phi_t*distmatT) * exp(- (phi_s*distmatS) / sqrt(1 + phi_t*distmatT))
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
                 epsilon = epsilon))
  
  elpd <- CV_sptvsampler(y = y, X = X, X_tilde = X_tilde,
                         N.samp = MC.samp,
                         V_z_full = V_z, L_z_full = L_z,
                         family = family, n_binom = n_binom,
                         mod_params = mod_params, CV_K = CV_K)
  
  return(list(beta = gamma[(n+1):(n+p), ],
              z = gamma[(n+p+1):(n+p+n*r), ],
              xi = gamma[1:n, ],
              elpd = elpd))
}

CV_sptvsampler <- function(y, X, X_tilde,
                           N.samp,
                           V_z_full, L_z_full = NULL,
                           family, n_binom = NULL,
                           mod_params, CV_K = 10){
  
  n <- length(y)
  partition_list <- id_partition(n, CV_K, random = FALSE)
  
  CV_samps <- lapply(1:length(partition_list), function(x)
    elpd_sptvGCM(y_train = y[-partition_list[[x]]],
                 X_train = X[-partition_list[[x]], ],
                 X_tilde_train = X_tilde[-partition_list[[x]], ],
                 y_pred = y[partition_list[[x]]],
                 X_pred = X[partition_list[[x]], ],
                 X_tilde_pred = X_tilde[partition_list[[x]], ],
                 J_tilde = V_z_full[- partition_list[[x]],
                                    partition_list[[x]]],
                 V_tilde = V_z_full[partition_list[[x]],
                                    partition_list[[x]]],
                 V_z_train = V_z_full[-partition_list[[x]],
                                      -partition_list[[x]]],
                 L_z_train = cholesky_CV(L_z_full, partition_list[[x]]),
                 N.samp = N.samp, mod_params = mod_params,
                 family = family,
                 n_binom_train = n_binom[-partition_list[[x]]],
                 n_binom_pred = n_binom[partition_list[[x]]]))
  
  elpd <- array(dim = n)
  for(k in 1:CV_K){
    elpd[partition_list[[k]]] = CV_samps[[k]]
  }
  return(elpd)
}

elpd_sptvGCM <- function(y_train, X_train, X_tilde_train,
                         y_pred, X_pred, X_tilde_pred, N.samp,
                         J_tilde, V_tilde, V_z_train, L_z_train = NULL,
                         family, mod_params,
                         n_binom_train = NULL, n_binom_pred = NULL){
  
  n <- length(y_train)
  p <- dim(X_train)[2]
  r <- dim(X_tilde_train)[2]
  
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  epsilon <- mod_params$epsilon
  
  # pre-allocation of useful Cholesky
  diag1 <- apply(X_tilde_train, 1, function(x){ 1/(2 + sum(x^2)) })
  S_D_star_train <- crossprod(X_train, diag1 * X_train)
  diag(S_D_star_train) <- diag(S_D_star_train) + 1
  chol_p_train <- chol(S_D_star_train)
  chol_nr_train <- lapply(1:n, function(x){
    chol(tcrossprod(X_tilde_train[x, ]) + diag(r))})
  
  # sample from the posterior
  gamma_train <- sapply(1:N.samp, function(x) 
    sampler_sptv(y = y_train, X = X_train, X_tilde = X_tilde_train,
                 family = family, n_binom = n_binom_train, 
                 chol_p = chol_p_train, chol_nr = chol_nr_train, 
                 L_z = L_z_train, nu_xi = nu_xi, nu_beta = nu_beta, 
                 nu_z = nu_z, epsilon = epsilon))
  
  samp_beta <- gamma_train[(n+1):(n+p), ]
  samp_z <- gamma_train[(n+p+1):(n+p+n*r), ]
  
  nk <- dim(X_tilde_pred)[1]
  z_pred <- array(dim = c(nk*r, N.samp))
  for(i in 1:r){
    ids <- ((i-1)*n+1):(i*n)
    idsk <- ((i-1)*nk+1):(i*nk)
    z_pred[idsk, ] <- predict_z(z_post = samp_z[ids, ], 
                                J = J_tilde, cholV = L_z_train,
                                V_tilde = V_tilde, nu_z = nu_z)
  }
  
  X_tilde_z <- sptv_prod(X_tilde_pred, z_pred)
  mu <- exp(X_pred %*% samp_beta + X_tilde_z)
  if(family == "poisson"){
    elpd_mat <- dpois(y_pred, mu, log = FALSE)
  }else if(family == "binomial"){
    prob <- ilogit(mu)
    elpd_mat <- dbinom(y_pred, n_binom_pred, prob, log = FALSE)
  }
  
  elpd <- apply(elpd_mat, 1, function(x) mean(x, na.rm = FALSE))
  
  return(log(elpd))
}

sptv_prod <- function(X_tilde, z){
  n <- dim(X_tilde)[1]
  r <- dim(X_tilde)[2]
  if(is.null(dim(z))){
    if(!(length(z) == (n*r))) stop("Dimension mismatch in sptv_prod().")
    out <- sapply(1:n, function(x) sum(X_tilde[x, ] * z[((1:r)-1)*n + x]))
  }else{
    if(!(dim(z)[1] == (n*r))) stop("Dimension mismatch in sptv_prod().")
    out <- t(sapply(1:n, function(x) X_tilde[x, ] %*% z[((1:r)-1)*n + x, ]))
  }
  return(out)
}

# n1 <- 50
# r1 <- 3
# X1 <- cbind(rep(1, n1), sapply(1:(r1-1), function(x) rnorm(n1)))
# z1 <- sapply(1:10, function(x) rnorm(n1 * r1))
# res1 <- sptv_prod(X1, z1)
# G1 <- makeG(X1)
# res2 <- G1 %*% z1
# summary(as.numeric(res1 - res2))

