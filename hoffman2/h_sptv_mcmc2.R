rm(list = ls())

source("../src/runsrc.R")

n_h <- 100
n_rep <- 2
samplesize_seq <- 1:5*100
# samplesize_seq <- 50
n_train_seq <- rep(samplesize_seq, each = n_rep)
n_postsamp_mcmc <- 5000

simdat <- read.csv("../data/sim_sptv_1000.csv")

temp <- simdat[1:4, ]
simdat[1:4, ] <- simdat[n_h + 1:4, ]
simdat[n_h + 1:4, ] <- temp

mlpd_out <- data.frame(matrix(ncol = 3, nrow = length(n_train_seq)))
names(mlpd_out) <- c("n", "mlpd", "runtime")
mlpd_out[, "n"] <- n_train_seq

write.csv(mlpd_out, "output/mlpd_poisson_mcmc2.csv", row.names = F)

burnin_pc <- 0.5
n_thin <- 5

for(k in 1:length(n_train_seq)){
  
  n_train <- n_train_seq[k]
  
  # Test on rows 1:n_h and train on next n_train
  simdat_h <- simdat[1:n_h, ]
  simdat_t <- simdat[n_h + 1:n_train, ]
  
  y_h <- as.numeric(simdat_h$y)
  X_h <- as.matrix(simdat_h[, grep("x", names(simdat_h))])
  X_tilde_h <- X_h
  S_h <- as.matrix(simdat_h[, c("s1", "s2")])
  time_h <- as.numeric(simdat_h$time)
  
  y <- as.numeric(simdat_t$y)
  X <- as.matrix(simdat_t[, grep("x", names(simdat_t))])
  X_tilde <- X
  S <- as.matrix(simdat_t[, c("s1", "s2")])
  time <- as.numeric(simdat_t$time)
  
  cat("Running n =", n_train, "...")
  t_mcmc_start <- Sys.time()
  mod_out <- sptvGLM_adaMetropGibbs2(y = y, X = X, X_tilde = X, 
                                     S = S, time = time, 
                                     family = "poisson", 
                                     N.samp = n_postsamp_mcmc, 
                                     starting = list(phi_s = 1, phi_t = 1,
                                                     beta = c(0, 0)), 
                                     prior = list(phi_s_a = 0.5, phi_s_b = 10,
                                                  phi_t_a = 0.5, phi_t_b = 10,
                                                  nu_xi = 1, nu_beta = 3,
                                                  nu_z = 3, alpha_epsilon = 0.5),
                                     n.batch = 3, batch.length = 10,
                                     verbose = F)
  t_mcmc_end <- Sys.time()
  
  mlpd_out[k, "runtime"] <- as.numeric(difftime(t_mcmc_end, t_mcmc_start), 
                                       units = "secs")
  
  ids <- 1:n_postsamp_mcmc
  ids <- ids[-(1:(floor(burnin_pc * n_postsamp_mcmc)))]
  ids <- ids[c(rep(FALSE, n_thin - 1), TRUE)]
  
  post_beta <- mod_out$beta[, ids]
  post_z <- mod_out$z[, ids]
  post_phi_s <- mod_out$phi_s[ids]
  post_phi_t <- mod_out$phi_t[ids]
  
  bigD_S <- as.matrix(dist(rbind(S_h, S)))
  bigD_t <- as.matrix(dist(c(time_h, time)))
  pd_mat <- array(dim = c(n_h, length(ids)))
  nk <- dim(X_tilde_h)[1]
  r <- dim(X_tilde_h)[2]
  
  cat("calculating MLPD")
  for(i in 1:length(ids)){
    
    phi_s <- post_phi_s[i]
    phi_t <- post_phi_t[i]
    nu_z <- 3                       # WARNING: SOFT CODE
    
    bigV <- 1/(1 + phi_t*bigD_t) * 
      exp(- (phi_s*bigD_S) / sqrt(1 + phi_t*bigD_t))
    chol_train <- chol(bigV[n_h + 1:n_train, n_h + 1:n_train])
    
    z_pred <- array(dim = c(nk*r))
    for(j in 1:r){
      ids <- ((j-1)*n_train+1):(j*n_train)
      idsk <- ((j-1)*nk+1):(j*nk)
      z_pred[idsk] <- predict_z(z_post = mod_out$z[ids, i],
                                J = bigV[n_h + 1:n_train, 1:n_h],
                                cholV = chol_train,
                                V_tilde = bigV[1:n_h, 1:n_h], nu_z = nu_z)
    }
    
    z_pred <- as.numeric(z_pred)
    X_tilde_h_z <- sptv_prod(X_tilde_h, z_pred)
    mu <- exp(X_h %*% post_beta + X_tilde_h_z)
    
    # y_pred <- t(apply(mu, 1, function(x) rpois(1:length(x), x)))
    p_y_pred <- t(sapply(1:n_h, function(x) dpois(y_h[x], mu[x])))
    p_y_pred <- apply(p_y_pred, 1, mean)
    
    pd_mat[, i] <- p_y_pred
  }
  
  pd_mc <- apply(pd_mat, 2, mean)
  mlpd <- mean(log(pd_mc))
  
  mlpd_out[k, "mlpd"] <- mlpd
  cat(" =", mlpd, "\n")
  
  write.csv(mlpd_out, "output/mlpd_poisson_mcmc2.csv", row.names = F)
}

