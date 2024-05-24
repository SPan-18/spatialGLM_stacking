# simulates count data with a latent GP on [0,1]x[0,1]
# with beta0 + beta1 * x1 + ... + beta_p * xp + z, z ~ 0.4 * GP(0, R_phi) 
# where x_k is a sample from N(0,1)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  # if(require(Rfast)){
  #   D <- Rfast::cholesky(V, parallel = TRUE)
  # }else{
  #   D <- chol(V)
  # }
  D <- chol(V)
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

sim_dyncount <- function(n, Nt, beta, phi, phi_t){
  
  S <- data.frame(s1 = c(0,0,1,1,runif(n*Nt - 4, 0, 1)),
                  s2 = c(0,1,0,1,runif(n*Nt - 4, 0, 1)))
  
  # p <- length(beta[[1]])
  p <- length(beta) / Nt
  
  X <- cbind(rep(1, n*Nt), sapply(1:(p-1), function(x) rnorm(n*Nt)))
  time <- rep(1:Nt, each = n)
  D <- as.matrix(dist(S))
  Dt <- as.matrix(dist(time))
  V <- 0.4 * exp(- phi * D) * exp(- phi_t * Dt)
  # V <- diag(n*Nt)
  
  z <- rmvn(1, rep(0, n*Nt), V)
  # z <- rnorm(n*Nt)
  # Xbeta <- array(dim = n*Nt)
  # for(i in 1:Nt){
  #   ids <- ((i-1)*n+1):(i*n)
  #   Xbeta[ids] <- X[ids, ] %*% beta[[i]]
  # }
  bigX <- mkdynX(X, time)
  Xbeta <- bigX %*% beta
  mu <- Xbeta + as.numeric(z)
  # mu <- X %*% beta + z + rnorm(n)
  y <- array(dim = n*Nt)
  for(i in 1:(n*Nt)){
    y[i] <- rpois(1, exp(mu[i]))
    # y[i] <- rnorm(1, mu[i], 1)
  }
  dat <- cbind(S, time, X, y = y, z = z)
  # dat <- cbind(S, time, X, y = rnorm(n, mu, 1), z = z)
  names(dat) = c("s1", "s2", "time", paste("x", 0:(p-1), sep = ""), "y", "z")
  return(dat)
}

# simdat <- sim_dyncount(n = 100, Nt = 10,
#                        beta = c(0.5, 0.7, 0.9, 1.2, 1.7, 2.2, 3, 4, 5.5, 7, 
#                                 -1.5, 0.2, 1.3, -0.3, 0.01, 0.2, -0.5, 0.12, -0.6, 0.26),
#                        phi = 3.5, phi_t = 0.5)
# write.csv(simdat, "../data/sim_dyncount100.10.csv", row.names = FALSE)

sim_sptcount2 <- function(n, Nt, beta, phi, phi_t){
  
  S <- data.frame(s1 = c(0,0,1,1,runif(n*Nt - 4, 0, 1)),
                  s2 = c(0,1,0,1,runif(n*Nt - 4, 0, 1)))
  
  # p <- length(beta[[1]])
  p <- length(beta)
  
  X <- cbind(rep(1, n*Nt), sapply(1:(p-1), function(x) rnorm(n*Nt)))
  time <- rep(1:Nt, each = n)
  D <- as.matrix(dist(S))
  Dt <- as.matrix(dist(time))
  V <- 0.4 * exp(- phi * D) * exp(- phi_t * Dt)
  z <- array(dim = n*Nt*p)
  for(i in 1:p){
    z[((i-1)*n*Nt+1):(i*n*Nt)] <- rmvn(1, rep(0, n*Nt), V)
  }
  G <- makeG(X)
  mu <- X %*% beta + G %*% z
  y <- array(dim = n*Nt)
  for(i in 1:(n*Nt)){
    y[i] <- rpois(1, exp(mu[i]))
  }
  dat <- cbind(S, time, X, y = y, z = z)
  names(dat) = c("s1", "s2", "time", paste("x", 0:(p-1), sep = ""), "y", "z")
  return(dat)
}

# simdat <- sim_sptcount(n = 100, Nt = 3,
#                        beta = c(5, -0.5),
#                        phi = 3.5, phi_t = 0.5)
# write.csv(simdat, "../data/sim_sptcount100.3.csv", row.names = FALSE)


sim_sptcount <- function(n, Nt, beta, phi, phi_t){
  
  S <- data.frame(s1 = c(0,0,1,1,runif(n*Nt - 4, 0, 1)),
                  s2 = c(0,1,0,1,runif(n*Nt - 4, 0, 1)))
  
  # p <- length(beta[[1]])
  p <- length(beta)
  
  X <- cbind(rep(1, n*Nt), sapply(1:(p-1), function(x) rnorm(n*Nt)))
  time <- rep(1:Nt, each = n)
  D <- as.matrix(dist(S))
  Dt <- as.matrix(dist(time))
  V <- 0.4 * exp(- phi * D) * exp(- phi_t * Dt)
  z <- rmvn(1, rep(0, n*Nt), V)
  mu <- X %*% beta + z
  y <- array(dim = n*Nt)
  for(i in 1:(n*Nt)){
    y[i] <- rpois(1, exp(mu[i]))
  }
  dat <- cbind(S, time, X, y = y, z = z)
  names(dat) = c("s1", "s2", "time", paste("x", 0:(p-1), sep = ""), "y", "z")
  return(dat)
}

# simdat <- sim_sptcount(n = 100, Nt = 3,
#                        beta = c(5, -0.5),
#                        phi = 3.5, phi_t = 0.5)
# write.csv(simdat, "../data/sim_sptmpcount100.3.csv", row.names = FALSE)

sim_sptvcount <- function(n, Nt, beta, phi_s, phi_t, sz){
  
  S <- vector(mode = "list", length = Nt)
  for(i in 1:Nt){
    S[[i]] <- data.frame(s1 = c(0,0,1,1,runif(n - 4, 0, 1)),
                         s2 = c(0,1,0,1,runif(n - 4, 0, 1)))
  }
  S <- do.call("rbind", S)
  
  p <- length(beta)
  
  X <- cbind(rep(1, n*Nt), sapply(1:(p-1), function(x) rnorm(n*Nt)))
  time <- rep(1:Nt, each = n)
  D <- as.matrix(dist(S))
  Dt <- as.matrix(dist(time))
  V <- 1/(1 + phi_t*Dt) * exp(- (phi_s*D) / sqrt(1 + phi_t*Dt))
  # z <- sapply(1:p, function(x){sz[x] * rmvn(1, rep(0, n*Nt), V)})
  # mu <- array(dim = n*Nt)
  # for(i in 1:length(mu)){
  #   mu[i] <- sum(X[i, ] * (beta + z[i, ]))
  # }
  # mu <- as.numeric(X %*% beta) + Xz
  L_z <- chol(V)
  z <- rnorm(n = n*Nt*p, mean = 0, sd = 1)
  for(j in 1:p){
    ids <- (j-1)*n*Nt+1:(n*Nt)
    z[ids] <- as.numeric(crossprod(L_z, z[ids]))
    z[ids] <- sqrt(sz[j]) * z[ids]
  }
  G <- makeG(X)
  mu <- (X %*% beta) + (G %*% z)
  y <- array(dim = n*Nt)
  for(i in 1:(n*Nt)){
    y[i] <- rpois(1, exp(mu[i]))
  }
  zmat <- matrix(z, nrow = n*Nt, ncol = p, byrow = F)
  dat <- cbind(S, time, X, y, zmat)
  names(dat) = c("s1", "s2", "time", paste("x", 0:(p-1), sep = ""), "y", 
                 paste("z", 1:p, sep = ""))
  return(dat)
}

# set.seed(1729)
# simdat <- sim_sptvcount(n = 1000, Nt = 3,
#                         beta = c(5, -0.5),
#                         phi_s = 3, phi_t = 0.5,
#                         sz = c(0.25, 0.25))
# write.csv(simdat, "../data/sim_sptvcount1000.3.csv", row.names = FALSE)

sim_sptv_continuous <- function(n, beta, phi_s, phi_t, sz){
  
  S <- data.frame(s1 = c(0,0,1,1,runif(n - 4, 0, 1)),
                  s2 = c(0,1,0,1,runif(n - 4, 0, 1)))
  time <- runif(n)
  p <- length(beta)
  
  X <- cbind(rep(1, n), sapply(1:(p-1), function(x) rnorm(n)))
  D <- as.matrix(dist(S))
  Dt <- as.matrix(dist(time))
  
  z <- rnorm(n = n*p, mean = 0, sd = 1)
  for(i in 1:p){
    V <- 1/(1 + phi_t[i]*Dt) * exp(- (phi_s[i]*D) / sqrt(1 + phi_t[i]*Dt))
    L_z <- chol(V)
    ids <- (i-1)*n+1:n
    z[ids] <- as.numeric(crossprod(L_z, z[ids]))
    z[ids] <- sqrt(sz[i]) * z[ids]
  }
  
  G <- makeG(X)
  mu <- (X %*% beta) + (G %*% z)
  y <- array(dim = n)
  for(i in 1:(n)){
    y[i] <- rpois(1, exp(mu[i]))
  }
  zmat <- matrix(z, nrow = n, ncol = p, byrow = F)
  dat <- cbind(S, time, X, y, zmat)
  names(dat) = c("s1", "s2", "time", paste("x", 0:(p-1), sep = ""), "y", 
                 paste("z", 1:p, sep = ""))
  return(dat)
}

# set.seed(1729)
# simdat <- sim_sptv_continuous(n = 5000, beta = c(5, -0.5),
#                               phi_s = c(2, 4), phi_t = c(0.5, 1),
#                               sz = c(0.25, 0.5))
# write.csv(simdat, "../data/sim_sptv_5000.csv", row.names = FALSE)
