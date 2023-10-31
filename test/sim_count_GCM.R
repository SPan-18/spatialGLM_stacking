rm(list = ls())

source("../src/runsrc.R")

# simulate count data
set.seed(1729)
N <- 500
p <- 2
y <- array(dim = N)
beta0 <- c(1, -1)
x.1 <- rnorm(N)
X <- cbind(rep(1, N), x.1)
Xbeta0 <- X %*% beta0
mu0 <- exp(Xbeta0)
for(i in 1:N) y[i] <- rpois(1, mu0[i])

# posterior of beta with CM prior; mu_beta = 0, V_beta = I_2
N.samp <- 1000
b <- 1
epsilon <- 0.75
s_beta_sq <- 5
s_xi_sq <- 5
post_beta_samp <- array(dim = c(N.samp, p))

H_gamma <- rbind(cbind(diag(N), X),
                 cbind(diag(N), matrix(0, N, p)),
                 cbind(matrix(0, p, N), diag(p)))
HtH <- t(H_gamma) %*% H_gamma
HtH.inv <- solve(HtH)

w_y <- rCM(n.samples = N.samp, mu = rep(0, N), V = diag(N),
           alpha = y + epsilon, kappa = rep(b, N), psi = "psi3")
w_beta <- rCM(n.samples = N.samp, mu = rep(0, p), V = diag(s_beta_sq, p),
              alpha = rep(10, p), kappa = rep(0.5, p), psi = "psi3")
w_xi <- rCM(n.samples = N.samp, mu = rep(0, N), V = diag(s_xi_sq, N),
            alpha = rep(0, N), kappa = rep(0.5, N), psi = "psi4")

w <- cbind(w_y, w_beta, w_xi)

# for(i in 1:N.samp){
#   w_y <- log(rgamma(N, y + epsilon, rep(b, N)))
#   w_beta <- (1/sqrt(rgamma(1,1,rgamma(1,100,1)))) * rnorm(p)
#   w_xi <- (1/sqrt(rgamma(1,1,rgamma(1,100,1)))) * rnorm(N)
#   w <- c(w_y, w_beta, w_xi)
#   
#   Htw <- t(H_gamma) %*% w
#   post_gamma <- HtH.inv %*% Htw
#   post_beta_samp[i, ] <- post_gamma[(N+1):(N+p)]
# }

Htw <- t(H_gamma) %*% t(w)
post_gamma_samp <- HtH.inv %*% Htw
post_gamma_samp <- t(post_gamma_samp)
post_beta_samp <- post_gamma_samp[, -c(1:N)]

print(ci_beta(post_beta_samp))


