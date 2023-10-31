# Sample a GCM posterior for a GLMM
# y = c(1, 2, 3), X = cbind(rep(1,3), c(1,2,3)), 
# S = data.frame(x = c(0.1, 0.2, 0.3), y = c(0.5, 0.9, 0.2))
# phi = 3, nu_xi = 0, nu_beta = 0.1, nu_z = 0.1, alpha_epsilon = 0.5

sampler_gaussian <- function(n, p, y, X, XtXplusIchol, L_z,
                             nu_xi, nu_beta, nu_z, alpha_epsilon){
  if(nu_xi == 0){
    sigmasq_xi <- 1
    w_xi <- rnorm(n)
  } else{
    sigmasq_xi <- 1/rgamma(1, 0.5 * nu_xi, 0.5 * nu_xi)
    w_xi <- rnorm(n) * sqrt(nu_xi / sigmasq_xi)
  } 
  sigmasq_beta <- 1/rgamma(1, 0.5 * nu_beta, 0.5 * nu_beta)
  sigmasq_z <- 1/rgamma(1, 0.5 * nu_z, 0.5 * nu_z)
  
  w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                     alpha = y[x] + alpha_epsilon, 
                                     kappa = 1, psi = "psi3"))
  w_beta <- rnorm(p) * sqrt(sigmasq_beta)
  w_z <- rnorm(n) * sqrt(sigmasq_z)
  w_z <- as.numeric(crossprod(L_z, w_z))
  gamma <- projection(X = X, prechol = XtXplusIchol, 
                      wy = w_y, wxi = w_xi, wbeta = w_beta, wz = w_z)
  return(gamma)
}

fast_posterior_GCM <- function(y, X, S = NULL, 
                               ypred = NULL, Xpred = NULL,
                               J_tilde = NULL, V_tilde = NULL,
                               V_z_train, N.samp,
                               family = "poisson",
                               beta_prior = "gaussian",
                               spatial_prior = "gaussian",
                               mod_params,
                               pred = FALSE,
                               Rfast_parallel = FALSE){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  phi <- mod_params$phi
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  
  L_z <- Rfast::cholesky(V_z_train, parallel = Rfast_parallel)
  
  XtXplusI <- crossprod(X) / 3 + diag(p)
  XtXplusIchol <- chol(XtXplusI)
  
  
  gamma <- sapply(1:N.samp, function(x) sampler_gaussian(n = n, p = p, y = y, X = X, 
                                                         XtXplusIchol = XtXplusIchol, 
                                                         L_z = L_z,
                                                         nu_xi = nu_xi, 
                                                         nu_beta = nu_beta, 
                                                         nu_z = nu_z, 
                                                         alpha_epsilon = alpha_epsilon))
  
  samp_beta <- gamma[(n+1):(n+p), ]
  samp_z <- gamma[(n+p+1):(2*n+p), ]
  # return(list(beta = gamma[(n+1):(n+p), ],
  #             z = gamma[(n+p+1):(2*n+p), ]))
  if(pred){
    if(is.null(ypred) && is.null(Xpred) && is.null(J_tilde) && is.null(V_tilde)){
      stop("Supply information on left out y, X, J_phi, V_phi.")
    }else{
      # find elpd
      z_tilde <- predict_z(z_post = gamma[(n+p+1):(2*n+p), ], 
                           J = J_tilde, cholV = L_z,
                           V_tilde = V_tilde, nu_z = nu_z)
      mu <- exp(Xpred %*% samp_beta + z_tilde)
      elpd_mat <- dpois(ypred, mu, log = TRUE)
      elpd <- apply(elpd_mat, 1, mean)
    }
    return(list(beta = samp_beta, z = samp_z, elpd = elpd))
  }else{
    return(list(beta = samp_beta,
                z = samp_z))
  }
}

posterior_GCM <- function(y, X, S, N.samp, 
                          family = "poisson",
                          fixed_prior = "gaussian", 
                          spatial_prior = "gaussian", 
                          mod_params){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  phi <- mod_params$phi
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  # kappa_epsilon <- mod_params$kappa_epsilon = 0
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = N.samp,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  distmat <- as.matrix(dist(S))
  XtX <- crossprod(X, X)
  HtH <- rbind(cbind(2*diag(n), X, diag(n)),
               cbind(t(X), XtX + diag(p), t(X)),
               cbind(diag(n), X, 2*diag(n)))
  chol11 <- Sys.time()
  HtHchol <- chol(HtH)
  chol12 <- Sys.time()
  chol21 <- Sys.time()
  L_z <- chol(exp(- phi * distmat))
  chol22 <- Sys.time()
  cat("Cholesky runtimes: \n", 
      "HtH:", chol12-chol11, units(chol12-chol11), "\n",
      "R(\U03C6):", chol22-chol21, units(chol22-chol21), "\n")
  beta_samps <- array(dim = c(N.samp, p))
  z_samps <- array(dim = c(N.samp, n))
  
  for(i in 1: N.samp){
    pb$tick()
    if(nu_xi == 0){
      sigmasq_xi <- 1
    } else{
      sigmasq_xi <- 1/rgamma(1, 1, nu_xi)
    } 
    sigmasq_beta <- 1/rgamma(1, 1, nu_beta)
    sigmasq_z <- 1/rgamma(1, 1, nu_z)
    
    w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                       alpha = y[x] + alpha_epsilon, 
                                       kappa = 1, psi = "psi3"))
    w_xi <- rnorm(n = n, mean = 0, sd = sqrt(sigmasq_xi))
    w_beta <- rnorm(n = p, mean = 0, sd = sqrt(sigmasq_beta))
    w_z <- rnorm(n = n, mean = 0, sd = sqrt(sigmasq_z))
    w_z <- as.numeric(crossprod(L_z, w_z))
    
    Htw <- c(w_y + w_xi, 
             as.numeric(crossprod(X, w_y)) + w_beta,
             w_y + w_z)
    temp <- forwardsolve(t(HtHchol), Htw)
    gamma <- backsolve(HtHchol, temp)
    beta_samps[i, ] <- gamma[(n+1):(n+p)]
    z_samps[i, ] <- gamma[(n+p+1):(2*n+p)]
  }
  return(list(beta = beta_samps,
              z = z_samps))
}

# # TEST: simulation example
# source("distributions.R")
# source("functions.R")
# source("projection.R")
# source("predictive.R")
# 
# library(MASS)
# library(progress)
# library(dplyr)
# library(ggplot2)
# library(gridExtra)
# # #
# set.seed(1729)
# n.grid <- 35
# Grid <- expand.grid(x.easting = seq(0, 1, length.out = n.grid),
#                     x.northing = seq(0, 1, length.out = n.grid))
# n <- nrow(Grid)
# holdout <- 1:10
# 
# # Explanatory variables and coefficients
# x1 <- rep(1, n)
# x2 <- rnorm(n) %>% round(2)
# 
# # Spatial field
# distance <- as.matrix(dist(Grid))
# phi <- 3
# Vz <- exp(- phi * distance)
# omega <- MASS::mvrnorm(n     = 1,
#                        mu    = rep(0,n),
#                        Sigma = 0.4 * Vz)
# beta0 <- 3
# beta1 <- 0.5
# eta <- beta0*x1 + beta1*x2 + omega
# 
# d <- Grid %>%
#   mutate(Y_pois   = rpois(n, exp(eta)),
#          x1       = x1,
#          x2       = x2)
# X <- as.matrix(cbind(x1, x2))
# #
# mod_params1 <- list(phi = 3, nu_xi = 0, nu_beta = 1000,
#                     nu_z = 1000, alpha_epsilon = 0.75)
# post.samp <- 100
# 
# t3 <- Sys.time()
# m.out_fast <- fast_posterior_GCM(y = d$Y_pois[-holdout], X = X[-holdout, ],
#                                  ypred = d$Y_pois[holdout],
#                                  Xpred = X[holdout, ],
#                                  J_tilde = Vz[-holdout, holdout],
#                                  V_tilde = Vz[holdout, holdout],
#                                  pred = TRUE,
#                                  V_z_train = Vz[-holdout, -holdout], 
#                                  N.samp = post.samp, mod_params = mod_params1)
# t4 <- Sys.time()
# print(t4 - t3)
# post_beta_fast <- t(m.out_fast$beta)
# print(ci_beta(post_beta_fast))
# post_z_fast <- m.out_fast$z
# postmean_z_fast <- apply(post_z_fast, 1, mean)
# 
# plot_simz <- ggplot(d %>% mutate(z=omega), aes(x = x.easting, y = x.northing)) +
#   geom_raster(aes(fill = z)) +
#   scale_fill_viridis_c(option = "inferno", direction = -1) +
#   # scale_fill_gradient2(midpoint = 2) +
#   xlab("Easting") +
#   ylab("Northing") +
#   theme_bw() +
#   theme(axis.line = element_line(color='black'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         aspect.ratio = 1)
# 
# plot_postz_fast <- ggplot(d[-holdout, ] %>% mutate(post_z=postmean_z_fast), aes(x = x.easting, y = x.northing)) +
#   geom_raster(aes(fill = post_z)) +
#   scale_fill_viridis_c(option = "inferno", direction = -1) +
#   # scale_fill_gradient2(midpoint = 2) +
#   xlab("Easting") +
#   ylab("Northing") +
#   theme_bw() +
#   labs(caption = paste("(based on", post.samp, "posterior samples)")) +
#   theme(axis.line = element_line(color='black'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         aspect.ratio = 1)
# 
# grid.arrange(plot_simz, plot_postz_fast, ncol=2)
