# Sample a GCM posterior for a GLMM
# y = c(1, 2, 3), X = cbind(rep(1,3), c(1,2,3)), 
# S = data.frame(x = c(0.1, 0.2, 0.3), y = c(0.5, 0.9, 0.2))
# phi = 3, nu_xi = 0, nu_beta = 0.1, nu_z = 0.1, alpha_epsilon = 0.5

sampler_GCM <- function(n, p, y, X, family,
                        n_binom = NULL,
                        XtXplusIchol, L_z,
                        nu_xi, nu_beta, nu_z, alpha_epsilon){
  # if(nu_xi == 0){
  #   sigmasq_xi <- 1
  #   w_xi <- rnorm(n)
  # } else{
  #   sigmasq_xi <- 1/rgamma(1, 0.5 * nu_xi, 0.5 * nu_xi)
  #   w_xi <- rnorm(n) * sqrt(nu_xi / sigmasq_xi)
  # } 
  w_xi <- rnorm(n, mean = 0, sd = nu_xi)
  sigmasq_beta <- 1/rgamma(1, 0.5 * nu_beta, 0.5 * nu_beta)
  sigmasq_z <- 1/rgamma(1, 0.5 * nu_z, 0.5 * nu_z)
  
  if(family == "poisson"){
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
  
  w_beta <- rnorm(p) * sqrt(sigmasq_beta)
  w_z <- rnorm(n) * sqrt(sigmasq_z)
  w_z <- as.numeric(crossprod(L_z, w_z))
  gamma <- projection(X = X, prechol = XtXplusIchol, 
                      wy = w_y, wxi = w_xi, wbeta = w_beta, wz = w_z)
  return(gamma)
}

elpd_GCM <- function(y_train, X_train, y_pred, X_pred, N.samp,
                     J_tilde, V_tilde, V_z_train, 
                     L_z_train = NULL,
                     family,
                     n_binom_train = NULL,
                     n_binom_pred = NULL,
                     beta_prior = "gaussian",
                     spatial_prior = "gaussian",
                     mod_params,
                     Rfastparallel = FALSE){
  
  n <- length(y_train)
  p <- dim(X_train)[2]
  
  # phi <- mod_params$phi
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  
  # attack here
  if(is.null(L_z_train)) L_z_train <- Rfast::cholesky(V_z_train, parallel = Rfastparallel)
  XtXplusI <- crossprod(X_train) / 3 + diag(p)
  XtXplusIchol <- chol(XtXplusI)
  
  gamma_train <- sapply(1:N.samp, function(x)  
    sampler_GCM(n = n, p = p, y = y_train, X = X_train, 
                family = family,
                n_binom = n_binom_train,
                XtXplusIchol = XtXplusIchol, 
                L_z = L_z_train,
                nu_xi = nu_xi, 
                nu_beta = nu_beta, 
                nu_z = nu_z, 
                alpha_epsilon = alpha_epsilon))
  
  # samp_xi <- gamma_train[1:n, ]
  samp_beta <- gamma_train[(n+1):(n+p), ]
  samp_z <- gamma_train[(n+p+1):(2*n+p), ]
  
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

slow_posterior_GCM <- function(y, X, S, N.samp, 
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
