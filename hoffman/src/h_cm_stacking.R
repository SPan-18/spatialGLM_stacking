# K-fold CV posterior sampling

CV_posterior_sampler <- function(y, X, N.samp,
                                 V_z_full, L_z_full = NULL,
                                 family,
                                 n_binom = NULL,
                                 beta_prior = "gaussian",
                                 spatial_prior = "gaussian",
                                 mod_params,
                                 CV_K = 10,
                                 Rfastparallel = FALSE){
  
  n <- length(y)
  partition_list <- id_partition(n, CV_K, random = FALSE)
  # phi <- mod_params$phi
  # V_z_full <- exp(- phi * distmat)
  # L_z_full <- Rfast::cholesky(V_z_full, parallel = Rfast_parallel)
  
  # ncores <- detectCores()
  CV_samps <- lapply(1:length(partition_list), function(x)
    elpd_GCM(y_train = y[-partition_list[[x]]],
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
             spatial_prior = spatial_prior, Rfastparallel = Rfastparallel))
  
  elpd <- array(dim = n)
  for(k in 1:CV_K){
    elpd[partition_list[[k]]] = CV_samps[[k]]
  }
  return(elpd)
}

# sample posterior and then find elpd

posterior_and_elpd <- function(y, X, N.samp, MC.samp,
                               distmat,
                               spCov = "matern",
                               n_binom = NULL,
                               family, 
                               beta_prior = "gaussian",
                               spatial_prior = "gaussian",
                               mod_params,
                               CV_K = 10,
                               Rfastparallel = FALSE){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  phi <- mod_params$phi
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  nu_matern <- mod_params$nu_matern
  
  if(spCov == "exponential"){
    V_z <- exp(- phi * distmat)
  }else if(spCov == "matern"){
    V_z <- matern(u = distmat, phi = 1/phi, kappa = nu_matern)
  }else{
    V_z <- geoR::cov.spatial(distmat, cov.model = spCov, 
                             cov.pars = c(1, 1/phi), kappa = nu_matern)
  }
  
  L_z <- Rfast::cholesky(V_z, parallel = Rfastparallel)
  XtXplusI <- crossprod(X) / 3 + diag(p)
  XtXplusIchol <- chol(XtXplusI)
  
  gamma <- sapply(1:N.samp, function(x) sampler_GCM(n = n, p = p, y = y, X = X,
                                                    family = family,
                                                    n_binom = n_binom,
                                                    XtXplusIchol = XtXplusIchol,
                                                    L_z = L_z,
                                                    nu_xi = nu_xi, 
                                                    nu_beta = nu_beta, 
                                                    nu_z = nu_z, 
                                                    alpha_epsilon = alpha_epsilon))
  
  # K-fold CV to find elpd
  elpd <- CV_posterior_sampler(y = y, X = X, N.samp = MC.samp,
                               V_z_full = V_z, L_z_full = L_z,
                               family = family,
                               n_binom = n_binom,
                               beta_prior = beta_prior,
                               spatial_prior = spatial_prior,
                               mod_params = mod_params,
                               CV_K = CV_K,
                               Rfastparallel = Rfastparallel)
  
  return(list(beta = gamma[(n+1):(n+p), ],
              z = gamma[(n+p+1):(2*n+p), ], 
              elpd = elpd))
}

# stacking using K-fold cross-validation

spGLM_stack <- function(y, X, S, N.samp, MC.samp = 200,
                         family = "poisson",
                         n_binom = NULL,
                         spCov = "matern",
                         beta_prior = "gaussian",
                         spatial_prior = "gaussian",
                         mod_params_list,
                         CV_fold = 10,
                         Rfastparallel = FALSE,
                         verbose = TRUE,
                         print_stackweights = TRUE){
  
  # distmat <- as.matrix(dist(S))
  distmat <- Rfast::Dist(S)
  
  n <- length(y)
  permut <- sample(1:n)
  y <- y[permut]
  if(!is.null(n_binom)){ n_binom <- n_binom[permut] }
  X <- X[permut, ]
  S <- S[permut, ]
  
  t_start <- Sys.time()
  samps <- lapply(1:length(mod_params_list), function(x)
    posterior_and_elpd(y = y, X = X,
                       distmat = distmat,
                       spCov = spCov,
                       N.samp = N.samp, MC.samp = MC.samp,
                       family = family,
                       n_binom = n_binom,
                       beta_prior = beta_prior,
                       spatial_prior = spatial_prior,
                       mod_params = mod_params_list[[x]],
                       CV_K = CV_fold, Rfastparallel = Rfastparallel))
  
  # samps <- vector(mode = "list", length = length(mod_params_list))
  # for(x in 1:length(mod_params_list)){
  #   samps[[x]] = posterior_and_elpd(y = y, X = X, 
  #                                   distmat = distmat,
  #                                   N.samp = N.samp, MC.samp = MC.samp,
  #                                   family = family,
  #                                   beta_prior = beta_prior,
  #                                   spatial_prior = spatial_prior,
  #                                   mod_params = mod_params_list[[x]],
  #                                   CV_K = CV_fold, Rfastparallel = Rfastparallel)
  # }
  
  elpd_mat <- do.call(cbind, lapply(samps, function(x) x$elpd))
  w_hat <- loo::stacking_weights(elpd_mat)
  t_end <- Sys.time()
  runtime <- difftime(t_end, t_start)
  if(verbose) cat("\nRUNTIME:", round(runtime, 2), units(runtime), ".\n\n")
  
  if(verbose){
    stack_out <- as.matrix(do.call(rbind, lapply(mod_params_list, unlist)))
    stack_out <- cbind(stack_out, round(as.numeric(w_hat), 2))
    colnames(stack_out) = c("phi", "smooth", "epsilon", "nu.xi", "nu.beta", "nu.z", "weight")
    rownames(stack_out) = paste("Model", 1:nrow(stack_out))
    if(print_stackweights){
      cat("MODEL WEIGHTS:\n")
      print(knitr::kable(stack_out))
    } 
  }
  
  for(i in 1:length(samps)){
    samps[[i]]$z <- samps[[i]]$z[order(permut), ]
    samps[[i]]$elpd <- samps[[i]]$elpd[order(permut)]
  }
  
  return(list(models = samps, weights = w_hat))
}

spGLM_onlystack <- function(y, X, S, MC.samp = 200,
                            family = "poisson",
                            n_binom = NULL,
                            spCov = "matern",
                            beta_prior = "gaussian",
                            spatial_prior = "gaussian",
                            mod_params_list,
                            CV_fold = 10,
                            Rfastparallel = FALSE,
                            verbose = TRUE,
                            print_stackweights = TRUE){
  
  distmat <- Rfast::Dist(S)
  n <- length(y)
  permut <- sample(1:n)
  y <- y[permut]
  if(!is.null(n_binom)){ n_binom <- n_binom[permut] }
  X <- X[permut, ]
  S <- S[permut, ]
  
  t_start <- Sys.time()
  elpd_list <- lapply(1:length(mod_params_list), function(x)
    CV_elpd(y = y, X = X,
            distmat = distmat,
            spCov = spCov,
            MC.samp = MC.samp,
            family = family,
            n_binom = n_binom,
            beta_prior = beta_prior,
            spatial_prior = spatial_prior,
            mod_params = mod_params_list[[x]],
            CV_K = CV_fold, Rfastparallel = Rfastparallel))
  
  elpd_mat <- do.call(cbind, elpd_list)
  w_hat <- loo::stacking_weights(elpd_mat)
  t_end <- Sys.time()
  runtime <- difftime(t_end, t_start)
  if(verbose) cat("\nRUNTIME:", round(runtime, 2), units(runtime), ".\n\n")
  
  if(verbose){
    stack_out <- as.matrix(do.call(rbind, lapply(mod_params_list, unlist)))
    stack_out <- cbind(stack_out, round(as.numeric(w_hat), 2))
    colnames(stack_out) = c("phi", "smooth", "epsilon", "nu.xi", "nu.beta", "nu.z", "weight")
    rownames(stack_out) = paste("Model", 1:nrow(stack_out))
    if(print_stackweights){
      cat("MODEL WEIGHTS:\n")
      print(knitr::kable(stack_out))
    } 
  }
  
  return(list(elpd = elpd_mat[order(permut), ],
              weights = w_hat))
}

CV_elpd <- function(y, X, MC.samp,
                    distmat,
                    spCov = "matern",
                    n_binom = NULL,
                    family, 
                    beta_prior = "gaussian",
                    spatial_prior = "gaussian",
                    mod_params,
                    CV_K = 10,
                    Rfastparallel = FALSE){
  
  phi <- mod_params$phi
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  nu_matern <- mod_params$nu_matern
  
  if(spCov == "exponential"){
    V_z <- exp(- phi * distmat)
  }else if(spCov == "matern"){
    V_z <- matern(u = distmat, phi = 1/phi, kappa = nu_matern)
  }else{
    V_z <- geoR::cov.spatial(distmat, cov.model = spCov, 
                             cov.pars = c(1, 1/phi), kappa = nu_matern)
  }
  L_z <- Rfast::cholesky(V_z, parallel = Rfastparallel)
  
  elpd <- CV_posterior_sampler(y = y, X = X, N.samp = MC.samp,
                               V_z_full = V_z, L_z_full = L_z,
                               family = family,
                               n_binom = n_binom,
                               beta_prior = beta_prior,
                               spatial_prior = spatial_prior,
                               mod_params = mod_params,
                               CV_K = CV_K,
                               Rfastparallel = Rfastparallel)
  return(elpd)
}

postrunsampler <- function(out, N.samp){
  stack_weights <- as.numeric(out$weights)
  n_post <- dim(out$models[[1]]$beta)[2]
  p.obs <- dim(out$models[[1]]$beta)[1]
  n.obs <- dim(out$models[[1]]$z)[1]
  ids <- sample(1:n_post, size = N.samp, replace = TRUE)
  post_samples <- sapply(1:N.samp, function(x){
    model_id <- sample(1:length(out$models), 1, prob = stack_weights)
    return(c(out$models[[model_id]]$beta[, ids[x]], 
             out$models[[model_id]]$z[, ids[x]]))
  })
  return(list(beta = post_samples[1:p.obs, ],
              z = post_samples[(p.obs+1):(n.obs+p.obs), ]))
}

# TEST
