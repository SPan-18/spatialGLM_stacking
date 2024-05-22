rm(list = ls())

source("../src/runsrc.R")

n_h <- 10
n_rep <- 20
samplesize_seq <- 1:5*100
# samplesize_seq <- 50
n_train_seq <- rep(samplesize_seq, each = n_rep)
# n_train_seq <- c(50, 50)
n_postsamp_stack <- 500

simdat <- read.csv("../data/sim_sptv_1000.csv")

mlpd_out <- data.frame(matrix(ncol = 3, nrow = length(n_train_seq)))
names(mlpd_out) <- c("n", "mlpd", "runtime")
mlpd_out[, "n"] <- n_train_seq

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
  
  mod_list <- list(G_phi_s = c(2, 3, 10), G_phi_t = c(0.25, 0.9, 1.5),
                   G_epsilon = c(0.25, 0.5), G_nu_xi = c(1),
                   G_nu_beta = c(3), G_nu_z = c(3))
  
  mod_list <- create_candidate_models(mod_list)
  
  t_stack_start <- Sys.time()
  m_out <- sptvGLM_stack(y = y, X = X, X_tilde = X, S = S, time = time,
                         N.samp = n_postsamp_stack, MC.samp = 200,
                         family = "poisson", mc.cores = 6, solver = "MOSEK",
                         mod_params_list = mod_list)
  t_stack_end <- Sys.time()
  
  # n_postsamp_mcmc <- (n_postsamp_stack * n_thin) / (1 - burnin_pc)
  # t_mcmc_start <- Sys.time()
  # mod_out <- sptvGLM_adaMetropGibbs(y = y, X = X, X_tilde = X, 
  #                                   S = S, time = time, 
  #                                   family = "poisson", 
  #                                   N.samp = n_postsamp_mcmc, 
  #                                   starting = list(phi_s = c(1, 1), 
  #                                                   phi_t = c(1, 1),
  #                                                   beta = c(0, 0)), 
  #                                   prior = list(phi_s_a = c(0.1, 0.1),
  #                                                phi_s_b = c(10, 10),
  #                                                phi_t_a = c(0.1, 0.1),
  #                                                phi_t_b = c(10, 10),
  #                                                nu_xi = 1, nu_beta = 3,
  #                                                nu_z = 3, alpha_epsilon = 0.5),
  #                                   n.batch = 3, batch.length = 10)
  # t_mcmc_end <- Sys.time()
  
  mlpd_out[k, "runtime"] <- as.numeric(difftime(t_stack_end, t_stack_start), units = "secs")
  
  bigD_S <- as.matrix(dist(rbind(S_h, S)))
  bigD_t <- as.matrix(dist(c(time_h, time)))
  L <- length(m_out$models)
  lpd_mat_stack <- array(dim = c(n_h, L))
  
  for(i in 1:L){
    cat("Stack lpd: n = ", n_train, " Model", i, "\n")
    phi_s <- mod_list[[i]]$phi_s
    phi_t <- mod_list[[i]]$phi_t
    nu_z <- mod_list[[i]]$nu_z
    
    bigV <- V_z <- 1/(1 + phi_t*bigD_t) * exp(- (phi_s*bigD_S) / sqrt(1 + phi_t*bigD_t))
    chol_train <- chol(bigV[n_h + 1:n_train, n_h + 1:n_train])
    
    nk <- dim(X_tilde_h)[1]
    r <- dim(X_tilde_h)[2]
    z_pred <- array(dim = c(nk*r, n_postsamp_stack))
    for(j in 1:r){
      ids <- ((j-1)*n_train+1):(j*n_train)
      idsk <- ((j-1)*nk+1):(j*nk)
      z_pred[idsk, ] <- predict_z(z_post = m_out$models[[i]]$z[ids, ], 
                                  J = bigV[n_h + 1:n_train, 1:n_h], 
                                  cholV = chol_train,
                                  V_tilde = bigV[1:n_h, 1:n_h], nu_z = nu_z)
    }
    
    X_tilde_h_z <- sptv_prod(X_tilde_h, z_pred)
    mu <- exp(X_h %*% m_out$models[[i]]$beta + X_tilde_h_z)
    
    y_pred <- t(apply(mu, 1, function(x) rpois(1:length(x), x)))
    p_y_pred <- t(sapply(1:n_h, function(x) dpois(y_pred[x, ], y_h[x])))
    p_y_pred <- apply(p_y_pred, 1, mean)
    lpd_mat_stack[, i] <- p_y_pred
  }
  
  w_opt <- m_out$weights
  mlpd <- mean(log(lpd_mat_stack %*% w_opt))
  
  mlpd_out[k, "mlpd"] <- mlpd
}

write.csv(mlpd_out, "output/mlpd_poisson_stack.csv", row.names = F)
