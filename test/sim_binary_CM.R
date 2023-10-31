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
p0 <- logit.inv(Xbeta0)
for(i in 1:N) y[i] <- rbinom(1, 1, p0[i])

# posterior of beta with CM prior; mu_beta = 0, V_beta = I_2
N.samp <- 100
b <- 1
epsilon <- 0.001
alpha_beta <- 0.1
kappa_beta <- 0.5
alpha_prior <- rep(alpha_beta, 2)
kappa_prior <- rep(kappa_beta, 2)
alpha_post <- c(0.5*y + 0.5*epsilon,
                0.5*y + 0.5*epsilon,
                alpha_prior)
kappa_post <- c(rep(b, N), rep(b, N), kappa_prior)
l_b <- 100

prior_beta_samp <- rCM(n.samples = N.samp, mu = c(0, 0), V = diag(rep(l_b, 2)), 
                       alpha = alpha_prior, kappa = kappa_prior, psi = "psi2")
prior_beta_dens <- kde2d(prior_beta_samp[,1], prior_beta_samp[,2], n = 50)
w_beta <- rCM(n.samples = N.samp, mu = rep(0, (2*N) + 2), V = diag((2*N) + 2),
              alpha = alpha_post, kappa = kappa_post, psi = "psi2")
H_beta <- rbind(X, X, diag(rep(1/l_b, 2)))
HtH <- t(H_beta) %*% H_beta
HtH.inv <- solve(HtH)
post_beta_samp <- HtH.inv %*% t(H_beta) %*% t(w_beta)
post_beta_samp <- t(post_beta_samp)
post_beta_dens <- kde2d(post_beta_samp[,1], post_beta_samp[,2], n = 50)
x_min = min(prior_beta_samp[, 1], post_beta_samp[, 1], beta0[1])
x_max = max(prior_beta_samp[, 1], post_beta_samp[, 1], beta0[1]+1)
y_min = min(prior_beta_samp[, 2], post_beta_samp[, 2], beta0[2])
y_max = max(prior_beta_samp[, 2], post_beta_samp[, 2], beta0[2])

contour(prior_beta_dens, lwd = 0.5, 
        # xlim = c(x_min, x_max), ylim = c(y_min, y_max),
        xlim = c(-600, 100), ylim = c(-500, 500),
        col = "lightblue4",
        main = bquote(alpha[beta] == .(alpha_beta)~","~
                        kappa[beta] == .(kappa_beta)))
contour(post_beta_dens, lwd = 0.8, add = TRUE, col = "red")
points(beta0[1], beta0[2], pch = 8, col = "darkgreen")

print(ci_beta(post_beta_samp))

print(HtH.inv %*% t(rbind(X, X)) %*% kbar_binary(c(0.5*y + 0.5*epsilon,
                                             0.5*y + 0.5*epsilon),
                                           c(rep(b, N), rep(b, N))))
