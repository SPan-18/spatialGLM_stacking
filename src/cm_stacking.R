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
  post_samples = list(beta = gamma[(n+1):(n+p), ],
                      z = gamma[(n+p+1):(2*n+p), ])
  
  # K-fold CV to find elpd
  elpd <- CV_posterior_sampler(y = y, X = X, V_z_full = V_z, N.samp = MC.samp,
                               family = family,
                               beta_prior = beta_prior,
                               spatial_prior = spatial_prior,
                               mod_params = mod_params,
                               CV_K = CV_K,
                               Rfastparallel = Rfastparallel)
  
  return(list(post_samples = post_samples, elpd = elpd))
}

# stacking using K-fold cross-validation

cm_stacking <- function(y, X, S, N.samp, MC.samp = 200,
                        family = "poisson",
                        beta_prior = "gaussian",
                        spatial_prior = "gaussian",
                        mod_params_list,
                        CV_fold = 10,
                        Rfast_parallel = FALSE){
  
  distmat <- as.matrix(dist(S))
  
  samps <- lapply(1:length(mod_params_list), function(x)
    posterior_and_elpd(y = y, X = X, distmat = distmat,
                       N.samp = N.samp, MC.samp = MC.samp,
                       family = family,
                       beta_prior = beta_prior,
                       spatial_prior = spatial_prior,
                       mod_params = mod_params_list[[x]],
                       CV_K = CV_fold,
                       Rfastparallel = Rfast_parallel))
  
  return(samps)
}

# TEST above function




