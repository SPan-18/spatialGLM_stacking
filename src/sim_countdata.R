# simulates count data with a latent GP on [0,1]x[0,1]
# with beta0 + beta1 * x1 + ... + beta_p * xp + z, z ~ 0.4 * GP(0, R_phi) 
# where x_k is a sample from N(0,1)

sim_count <- function(n, beta, phi){
  S <- data.frame(x = c(0,0,1,1,runif(n - 4, 0, 1)),
                  y = c(0,1,0,1,runif(n - 4, 0, 1)))
  p <- length(beta)
  X <- cbind(rep(1, n), sapply(1:(p-1), function(x) rnorm(n)))
  V <- as.matrix(dist(S))
  z <- MASS::mvrnorm(n = 1, mu = rep(0,n), 
                     Sigma = 0.4 * exp(- phi * V))
  mu <- X %*% beta + z
  dat <- cbind(S, X, dat = rpois(n, exp(mu)))
  names(dat) = c("easting", "northing", paste("x", 0:(p-1), sep = ""), "y")
  return(dat)
}

# TEST
# sim_count(10, c(3, 1, -0.4), 3)
