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

CV_posterior_sampler <- function(y, X, distmat, N.samp,
                                 family = "poisson",
                                 beta_prior = "gaussian",
                                 spatial_prior = "gaussian",
                                 mod_params,
                                 CV_K = 10,
                                 Rfast_parallel = FALSE){
  
  n <- length(y)
  partition_list <- id_partition(n, CV_K)
  phi <- mod_params$phi
  V_z_full <- exp(- phi * distmat)
  # L_z_full <- Rfast::cholesky(V_z_full, parallel = Rfast_parallel)
  
  CV_samps <- lapply(1:length(partition_list), function(x)
    fast_posterior_GCM(y = y[-partition_list[[x]]], 
                       X = X[-partition_list[[x]], ], 
                       ypred = y[partition_list[[x]]], 
                       Xpred = X[partition_list[[x]], ], 
                       J_tilde = V_z_full[- partition_list[[x]], 
                                          partition_list[[x]]],
                       V_tilde = V_z_full[partition_list[[x]], 
                                          partition_list[[x]]],
                       pred = TRUE,
                       V_z_train = V_z_full[-partition_list[[x]], 
                                            -partition_list[[x]]], 
                       N.samp = N.samp, 
                       mod_params = mod_params,
                       family = family,
                       beta_prior = beta_prior, 
                       spatial_prior = spatial_prior))
  elpd <- array(dim = n)
  for(k in 1:CV_K){
    elpd[partition_list[[k]]] = CV_samps[[k]]$elpd
  }
  return(list(samples = CV_samps, elpd = elpd))
}

# stacking using K-fold cross-validation

cm_stacking <- function(y, X, S, N.samp,
                        family = "poisson",
                        beta_prior = "gaussian",
                        spatial_prior = "gaussian",
                        mod_params_list,
                        CV_K = 10,
                        Rfast_parallel = FALSE){
  
  distmat <- as.matrix(dist(S))
  
  samps <- lapply(1:length(mod_params_list), function(x) 
    CV_posterior_sampler(y = y, X = X, distmat = distmat, N.samp = N.samp,
                         family = family,
                         beta_prior = beta_prior,
                         spatial_prior = spatial_prior,
                         mod_params = mod_params_list[[x]],
                         CV_K = 10,
                         Rfast_parallel = FALSE))
  return(samps)
}

# TEST above function




