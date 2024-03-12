rm(list = ls())

source("../src/runsrc.R")

n_h <- 100
n_train_seq <- c(100, 200, 300, 400, 500)
# n_train_seq <- c(50, 50)

simdat <- read.csv("../data/sim_count1000.csv")

mlpd_mat <- array(dim = c(length(n_train_seq), 2))
mlpd_mat[, 1] <- n_train_seq

bigD <- as.matrix(dist(rbind(S_h, S)))

for(k in 1:length(n_train_seq)){
  
  n_train <- n_train_seq[k]
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
  
  mod_list <- create_model_list(G_decay = c(3, 4, 10), 
                                G_smoothness = c(0.5, 1, 1.5),
                                G_epsilon = c(0.25, 0.5),
                                G_nuxi = 1,
                                G_nubeta = 2.1, G_nuz = 2.1)
  
  m_out <- spGLM_stack(y = y, X = X, S = S, N.samp = n_postsamp,
                       family = "poisson",
                       spCov = "matern",
                       solver = "MOSEK",
                       mc.cores = 6,
                       mod_params_list = mod_list)
  
  L <- length(m_out$models)
  lpd_mat <- array(dim = c(n_h, L))
  
  for(i in 1:L){
    cat("n = ", n_train, " Model", i, "\n")
    phi <- mod_list[[i]]$phi
    nu <- mod_list[[i]]$nu_matern
    nu_z <- mod_list[[i]]$nu_z
    
    bigV <- matern(u = bigD, phi = 1/phi, kappa = nu)
    chol_train <- chol(bigV[n_h + 1:n_train, n_h + 1:n_train])
    z_tilde <- predict_z(z_post = m_out$models[[i]]$z, 
                         J = bigV[n_h + 1:n_train, 1:n_h],
                         cholV = chol_train,
                         V_tilde = bigV[1:n_h, 1:n_h],
                         nu_z = nu_z)
    z_tilde[which(z_tilde < -500)] <- -500
    z_tilde[which(z_tilde > 500)] <- 500
    
    mu <- exp(X_h %*% m_out$models[[i]]$beta + z_tilde)
    y_pred <- t(apply(mu, 1, function(x) rpois(1:length(x), x)))
    p_y_pred <- t(sapply(1:n_h, function(x) dpois(y_pred[x, ], y_h[x])))
    p_y_pred <- apply(p_y_pred, 1, mean)
    lpd_mat[, i] <- p_y_pred
  }
  
  w_opt <- m_out$weights
  mlpd <- mean(log(lpd_mat %*% w_opt))
  
  mlpd_mat[k, 2] <- mlpd
}

write.csv(mlpd_mat, "mlpd_poisson_stack.csv", row.names = F)

