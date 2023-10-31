rm(list = ls())

source("../src/dist.R")
source("../src/functions.R")

# simulate count data
set.seed(1729)
N <- 100
y <- array(dim = N)
beta0 <- c(1, -1)
x.1 <- rnorm(N)
X <- cbind(rep(1, N), x.1)
Xbeta0 <- X %*% beta0
mu0 <- exp(Xbeta0)
for(i in 1:N) y[i] <- rpois(1, mu0[i])

# posterior of beta with CM prior; mu_beta = 0, V_beta = I_2
N.samp <- 100
b <- 1

prior.list <- list(alpha_beta = 100,
                   kappa_beta = 1,
                   l_b = 1,
                   mu = rep(0, 2),
                   epsilon = 0.05)

post_beta_samp <- posterior_psi3(y = y, X = X, b = b,
                                 n.samples = N.samp,
                                 prior = prior.list)

bseq <- c(0.5, 1, 1.5)
K = 10
partition_list <- id_partition(N, K)

elpd_hat <- array(dim = c(N, length(bseq)))
beta_k <- vector(mode = "list", length = K)

for(j in 1:length(bseq)){
  for(k in 1:K){
    beta_k[[k]][[j]] <- posterior_psi3(y = y[-partition_list[[k]]], 
                                       X = X[-partition_list[[k]], ], 
                                       b = bseq[j],
                                       n.samples = N.samp,
                                       prior = prior.list)
    # elpd
    p_hat <- array(dim = length(partition_list[[k]]))
    for(i in 1:length(partition_list[[k]])){
      mu_i <- exp(X[partition_list[[k]][i], ] %*% t(beta_k[[k]][[j]]))
      p_hat[i] <- mean(dpois(y[partition_list[[k]][i]], lambda = bseq[j] * mu_i))
    }
    elpd_hat[partition_list[[k]], j] <- log(p_hat)
  }
}

w_hat <- loo::stacking_weights(elpd_hat)
print(w_hat)

# J = 1
# s_J <- which(sapply(partition_list, function(x) (J %in% x)))
# mu_J <- exp(X[J, ] %*% t(beta_k[[s_J]][[1]]))
# ypred_J <- array(dim = length(mu_J))
# for(i in 1:length(mu_J)){
#   ypred_J[i] <- rpois(1, lambda = mu_J[i])
# }
# plot(density(ypred_J))


# post_beta_dens <- kde2d(post_beta_samp[,1], post_beta_samp[,2], n = 50)
# contour(post_beta_dens, lwd = 0.8,
#         xlim = c(-6, 5), ylim = c(-8, 5), col = "red")
# points(beta0[1], beta0[2], pch = 8, col = "darkgreen")