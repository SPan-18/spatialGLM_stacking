rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_count1000.csv")

n_h <- 100
n_train_seq <- c(100, 200, 300, 400, 500)
# n_train_seq <- c(50, 100, 200)

mlpd_mat <- array(dim = c(length(n_train_seq), 2))
mlpd_mat[, 1] <- n_train_seq

for(k in 1:length(n_train_seq)){
  
  n_train <- n_train_seq[k]
  cat("n = ", n_train, "\n")
  # Test on rows 1:n_h and train on next n_train
  simdat_h <- simdat[1:n_h, ]
  simdat_t <- simdat[n_h + 1:n_train, ]
  
  y_h <- as.numeric(simdat_h$y)
  X_h <- as.matrix(simdat_h[, grep("x", names(simdat_h))])
  S_h <- as.matrix(simdat_h[, c("s1", "s2")])
  
  y <- as.numeric(simdat_t$y)
  X <- as.matrix(simdat_t[, grep("x", names(simdat_t))])
  S <- as.matrix(simdat_t[, c("s1", "s2")])
  
  n_postsamp <- 1000
  
  mod_out <- spGCM_adaMetropGibbs(y = y, X = X, S = S, 
                                  family = "poisson", 
                                  N.samp = n_postsamp,
                                  spCov = "matern",
                                  starting = list(phi = 3, nu = 1,
                                                  beta = c(0, 0)),
                                  prior = list(phi = c(0.5, 10),
                                               nu = c(0.1, 2),
                                               nu_xi = 1, 
                                               nu_beta = 2.1,
                                               nu_z = 3,
                                               alpha_epsilon = 0.5))
  
  ids <- 1:n_postsamp
  ids <- ids[-(1:(floor(0.1 * n_postsamp))+1)]
  ids <- ids[c(rep(FALSE, 8), TRUE)]
  # ids <- sample(ids, size = 1000, replace = FALSE)
  
  post_beta <- mod_out$beta[, ids]
  post_z <- mod_out$z[, ids]
  post_phi <- mod_out$phi[ids]
  post_nu <- mod_out$nu[ids]
  
  bigD <- as.matrix(dist(rbind(S_h, S)))
  pd_mat <- array(dim = c(n_h, length(ids)))
  
  for(i in 1:length(ids)){
    
    phi <- post_phi[i]
    nu <- post_nu[i]
    nu_z <- 3
    
    bigV <- matern(u = bigD, phi = 1/phi, kappa = nu)
    chol_train <- chol(bigV[n_h + 1:n_train, n_h + 1:n_train])
    z_tilde <- predict_z(z_post = post_z[, i], 
                         J = bigV[n_h + 1:n_train, 1:n_h],
                         cholV = chol_train,
                         V_tilde = bigV[1:n_h, 1:n_h],
                         nu_z = nu_z)
    
    mu <- exp(X_h %*% post_beta[, i] + z_tilde)
    y_pred <- rpois(1:length(mu), mu)
    p_y_pred <- sapply(1:n_h, function(x) dpois(y_pred[x], y_h[x]))
    pd_mat[, i] <- p_y_pred
  }
  
  pd_mc <- apply(pd_mat, 1, mean)
  mlpd <- mean(log(pd_mc))
  cat("n = ", n_train, " mlpd = ", mlpd, "\n")
  
  mlpd_mat[k, 2] <- mlpd
}

write.csv(mlpd_mat, "mlpd_poisson_mcmc.csv", row.names = F)
