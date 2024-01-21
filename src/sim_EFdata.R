# simulates count data with a latent GP on [0,1]x[0,1]
# with beta0 + beta1 * x1 + ... + beta_p * xp + z, z ~ 0.4 * GP(0, R_phi) 
# where x_k is a sample from N(0,1)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  if(require(Rfast)){
    D <- Rfast::cholesky(V, parallel = TRUE)
  }else{
    D <- chol(V)
  }
  # D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p) %*% D + rep(mu,rep(n,p)))
}

sim_count <- function(n, beta, phi){
  S <- data.frame(s1 = c(0,0,1,1,runif(n - 4, 0, 1)),
                  s2 = c(0,1,0,1,runif(n - 4, 0, 1)))
  # S <- data.frame(s1 = runif(n, 0, 1),
  #                 s2 = runif(n, 0, 1))
  p <- length(beta)
  X <- cbind(rep(1, n), sapply(1:(p-1), function(x) rnorm(n)))
  D <- as.matrix(dist(S))
  V <- 0.4 * exp(- phi * D)
  z <- rmvn(1, rep(0, n), V)
  # z <- MASS::mvrnorm(n = 1, mu = rep(0,n), 
  #                    Sigma = 0.4 * exp(- phi * V))
  mu <- X %*% beta + z
  # mu <- X %*% beta + z + rnorm(n)
  dat <- cbind(S, X, y = rpois(n, exp(mu)), z = z)
  names(dat) = c("s1", "s2", paste("x", 0:(p-1), sep = ""), "y", "z")
  return(dat)
}

# TEST
# simdat <- sim_count(1000, c(5, -0.5), 3.5)
# write.csv(simdat, "../data/sim_count1000v2.csv")

ilogit <- function(x){
  return(1.0 / (1.0 + exp(- x)))
}

sim_binom <- function(n, n_binom, beta, phi){
  S <- data.frame(s1 = c(0,0,1,1,runif(n - 4, 0, 1)),
                  s2 = c(0,1,0,1,runif(n - 4, 0, 1)))
  p <- length(beta)
  X <- cbind(rep(1, n), sapply(1:(p-1), function(x) rnorm(n)))
  V <- as.matrix(dist(S))
  z <- MASS::mvrnorm(n = 1, mu = rep(0,n), 
                     Sigma = 0.4 * exp(- phi * V))
  mu <- X %*% beta + z
  prob <- ilogit(mu)
  y2 <- rpois(n, n_binom)
  dat <- cbind(S, X, y1 = rbinom(n, y2, prob = prob), y2 = y2, z = z)
  names(dat) = c("s1", "s2", paste("x", 0:(p-1), sep = ""), "y1", "y2", "z")
  return(dat)
}

sim_binom_nonsp <- function(n, n_binom, beta){
  p <- length(beta)
  X <- cbind(rep(1, n), sapply(1:(p-1), function(x) rnorm(n)))
  mu <- X %*% beta
  prob <- ilogit(mu)
  dat <- cbind(X, y = rbinom(n, rep(n_binom, n), prob = prob))
  colnames(dat) = c(paste("x", 0:(p-1), sep = ""), "y")
  return(as.data.frame(dat))
}

# TEST
# simdat = sim_binom(n = 1000, n_binom = 20, beta = c(1, -0.5), phi = 3.5)
# write.csv(simdat, "../data/sim_binom1000.csv", row.names = FALSE)
