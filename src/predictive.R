# (MVT) predictive spatial process multivariate t

predict_z <- function(z_post, J, cholV, V_tilde, nu_z, 
                      Rfast_parallel = FALSE){
  m <- dim(V_tilde)[1]
  n.samp <- dim(z_post)[2]
  n <- dim(z_post)[1]
  z_tilde <- backsolve(cholV, z_post, transpose = TRUE)
  # browser()
  u <- rgamma(n.samp, 0.5 * (nu_z + n), 0.5)
  w <- sapply(1:n.samp, function(x) rnorm(m) * sqrt((nu_z + sum(z_tilde^2)) / u[x]))
  # w <- sapply(1:100, function(x) rnorm(m))
  temp <- backsolve(cholV, J, transpose = TRUE)
  z_tilde <- crossprod(temp, z_tilde)
  # return(z_tilde)
  temp <- crossprod(temp)
  temp <- as.matrix(V_tilde) - temp
  temp <- Rfast::cholesky(temp, parallel = Rfast_parallel)
  z_tilde <- z_tilde + crossprod(temp, w)
  return(z_tilde)
}

### Do for multiple z

# TEST predict_z function
# n <- 100
# n_post <- 100
# m <- 3
# S <- data.frame(x.easting = runif(n, 0, 1),
#                 x.northing = runif(n, 0, 1))
# dist_S <- as.matrix(dist(S))
# phi <- 3
# nu_z <- 100
# V <- exp(- phi * dist_S)
# Vchol <- Rfast::cholesky(V)
# u <- rgamma(n_post, 0.5 * nu_z, 0.5)
# w <- sapply(1:n_post, function(x) rnorm(n) * sqrt(nu_z / u[x]))
# # w <- rnorm(n)
# z <- crossprod(Vchol, w)
# ids <- 1:m
# z_obs <- z[- ids, ]
# Vchol_obs <- Rfast::cholesky(V[- ids, - ids])
# 
# z_pred <- predict_z(z_post = z_obs, J = V[- ids, ids], cholV = Vchol_obs,
#                     V_tilde = V[ids, ids], nu_z = nu_z)
# 
# apply(z[ids,], 1, mean)
# apply(z_pred, 1, mean)
# 
# apply(t(V[- ids, ids]) %*% solve(V[- ids, - ids]) %*% z_obs, 1, summary)
