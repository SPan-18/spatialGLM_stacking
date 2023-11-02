# G_phi <- c(2, 3)
# G_alphaepsilon <- c(0.5, 0.75)
# model_list <- create_model_list(G_phi, G_alphaepsilon)

create_model_list <- function(G_phi, G_epsilon, G_nuxi = 0,
                              G_nubeta = 100, G_nuz = 100){
  models <- expand.grid(G_phi, G_epsilon, G_nuxi, G_nubeta, G_nuz)
  names(models) <- c("phi", "alpha_epsilon", "nu_xi", "nu_beta", "nu_z")
  model_list <- lapply(1:nrow(models), function(x){
    return(list(phi = models[x, "phi"],
                alpha_epsilon = models[x, "alpha_epsilon"],
                nu_xi = models[x, "nu_xi"],
                nu_beta = models[x, "nu_beta"],
                nu_z = models[x, "nu_z"]))})
  return(model_list)
}

# K-fold CV posterior sampling

CV_posterior_sampler <- function(y, X, V_z_full, N.samp,
                                 family = "poisson",
                                 beta_prior = "gaussian",
                                 spatial_prior = "gaussian",
                                 mod_params,
                                 CV_K = 10,
                                 Rfastparallel = FALSE){
  
  n <- length(y)
  partition_list <- id_partition(n, CV_K)
  # phi <- mod_params$phi
  # V_z_full <- exp(- phi * distmat)
  # L_z_full <- Rfast::cholesky(V_z_full, parallel = Rfast_parallel)
  
  ncores <- detectCores()
  CV_samps <- mclapply(1:length(partition_list), function(x)
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
             N.samp = N.samp, 
             mod_params = mod_params,
             family = family,
             beta_prior = beta_prior, 
             spatial_prior = spatial_prior,
             Rfastparallel = Rfastparallel),
    mc.cores = ncores)
  
  elpd <- array(dim = n)
  for(k in 1:CV_K){
    elpd[partition_list[[k]]] = CV_samps[[k]]
  }
  return(elpd)
}

# sample posterior and then find elpd

posterior_and_elpd <- function(y, X, distmat, N.samp, MC.samp,
                               family = "poisson", 
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
  
  V_z <- exp(- phi * distmat)
  L_z <- Rfast::cholesky(V_z, parallel = Rfastparallel)
  XtXplusI <- crossprod(X) / 3 + diag(p)
  XtXplusIchol <- chol(XtXplusI)
  
  gamma <- sapply(1:N.samp, function(x) sampler_poisson(n = n, p = p, y = y, X = X, 
                                                        XtXplusIchol = XtXplusIchol, 
                                                        L_z = L_z,
                                                        nu_xi = nu_xi, 
                                                        nu_beta = nu_beta, 
                                                        nu_z = nu_z, 
                                                        alpha_epsilon = alpha_epsilon))
  
  # K-fold CV to find elpd
  elpd <- CV_posterior_sampler(y = y, X = X, V_z_full = V_z, N.samp = MC.samp,
                               family = family,
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

cm_stacking <- function(y, X, S, N.samp, MC.samp = 200,
                        family = "poisson",
                        beta_prior = "gaussian",
                        spatial_prior = "gaussian",
                        mod_params_list,
                        CV_fold = 10,
                        Rfast_parallel = FALSE,
                        verbose = TRUE,
                        print_stackweights = TRUE){
  
  distmat <- as.matrix(dist(S))
  
  t_start <- Sys.time()
  samps <- lapply(1:length(mod_params_list), function(x)
    posterior_and_elpd(y = y, X = X, distmat = distmat,
                       N.samp = N.samp, MC.samp = MC.samp,
                       family = family,
                       beta_prior = beta_prior,
                       spatial_prior = spatial_prior,
                       mod_params = mod_params_list[[x]],
                       CV_K = CV_fold,
                       Rfastparallel = Rfast_parallel))
  
  elpd_mat <- do.call(cbind, lapply(samps, function(x) x$elpd))
  w_hat <- loo::stacking_weights(elpd_mat)
  t_end <- Sys.time()
  runtime <- difftime(t_end, t_start)
  if(verbose) cat("\nRUNTIME:", round(runtime, 2), units(runtime), ".\n\n")
  
  stack_out <- as.matrix(do.call(rbind, lapply(mod_params_list, unlist)))
  stack_out <- cbind(stack_out, round(as.numeric(w_hat), 2))
  colnames(stack_out) = c("phi", "epsilon", "nu.xi", "nu.beta", "nu.z", "weight")
  rownames(stack_out) = paste("Model", 1:nrow(stack_out))
  if(print_stackweights){
    cat("MODEL WEIGHTS:\n")
    print(knitr::kable(stack_out))
  } 
  return(list(models = samps, weights = w_hat))
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

# TEST above function




