# logit_bounds <- function(x, a, b){
#   return(log(x - a) - log(b - x))
# }
# 
# logitinv_bounds <- function(x, a, b){
#   return(a + (b - a) * (1 / (1 + exp(- x))))
#   # return(b - (b - a) / (1 + exp(x)))
# }

# rmvn <- function(n, mu=0, V = matrix(1)){
#   p <- length(mu)
#   if(any(is.na(match(dim(V),p))))
#     stop("Dimension problem!")
#   D <- chol(V)
#   t(matrix(rnorm(n*p), ncol=p) %*% D + rep(mu,rep(n,p)))
# }

# spGCM_aMCMC <- function(y, X, S, N.samp, family, 
#                         # spCov,
#                         n_binom = NULL,
#                         marginalPrior = TRUE,
#                         starting,
#                         prior,
#                         tuning,
#                         n.batch,
#                         batch.length,
#                         n.report,
#                         Rfastparallel = FALSE){
#   
#   n <- dim(X)[1]
#   p <- dim(X)[2]
#   
#   # prior for spatial process parameters
#   phiUnifa <- prior$phi[1]
#   phiUnifb <- prior$phi[2]
#   nuUnifa <- prior$nu[1]
#   nuUnifb <- prior$nu[2]
#   
#   # prior for variance parameters
#   nu_xi <- prior$nu_xi
#   nu_beta <- prior$nu_beta
#   nu_z <- prior$nu_z
#   
#   # boundary adjustment parameter
#   alpha_epsilon <- prior$alpha_epsilon
#   
#   H <- rbind(cbind(diag(n), X, diag(n)),
#              cbind(diag(n), matrix(0, n, p), matrix(0, n, n)),
#              cbind(matrix(0, p, n), diag(p), matrix(0, p, n)),
#              cbind(matrix(0, n, n), matrix(0, n, p), diag(n)))
#   XtX <- crossprod(X, X)
#   HtH <- rbind(cbind(2*diag(n), X, diag(n)),
#                cbind(t(X), XtX + diag(p), t(X)),
#                cbind(diag(n), X, 2*diag(n)))
#   HtHchol <- Rfast::cholesky(HtH, parallel = Rfastparallel)
#   PH_sqrt <- backsolve(HtHchol, t(H), transpose = TRUE)
#   PH <- crossprod(PH_sqrt)
#   QQt <- diag(3*n + p) - PH
#   eig_QQt <- eigen(QQt)
#   Q <- eig_QQt$vectors[, 1:n]
#   
#   distmat <- Rfast::Dist(S)
#   
#   target <- function(theta){
#     xi <- theta[1:n]
#     beta <- theta[(n+1):(n+p)]
#     z <- theta[(n+p+1):(2*n+p)]
#     sigmabetasq <- exp(theta[2*n+p+1])
#     sigmazsq <- exp(theta[2*n+p+2])
#     phi <- logitinv_bounds(theta[2*n+p+3], phiUnifa, phiUnifb)
#     nu <- logitinv_bounds(theta[2*n+p+4], nuUnifa, nuUnifb)
#     # q <- theta[(2*n+p+5):(3*n+p+4)]
#     
#     Vz <- matern(u = distmat, phi = 1/phi, kappa = nu)
#     Lz <- Rfast::cholesky(Vz, parallel = Rfastparallel)
#     Qz <- sum(backsolve(Lz, z, transpose = TRUE)^2)
#     Qbeta <- sum(beta^2)
#     # Qq <- Q %*% q
#     # mu <- - 1.0 * Qq[1:n]
#     logLambda <- X %*% beta + z + xi #- mu
#     if(family == "poisson"){
#       Qxi <- sum(xi^2)
#       lLik <- sum(dpois(y, exp(logLambda), log = TRUE))
#     }
#     
#     logPost <- (
#       ##Priors
#       -(0.5 * nu_beta + 1) * log(sigmabetasq) - (0.5 * nu_beta)/sigmabetasq
#       -(0.5 * nu_z + 1) * log(sigmazsq) - (0.5 * nu_z)/sigmazsq
#       - 0.5 * n * log(sigmazsq) - 0.5 * sum(log(diag(Lz))) - 0.5 * Qz/sigmazsq
#       - 0.5 * p * log(sigmabetasq) - 0.5 * Qbeta/sigmabetasq
#       # + sum(alpha_epsilon * logLambda) 
#       - 0.5 * Qxi / (nu_xi)^2
#       ##Jacobians
#       +log(sigmabetasq) + log(sigmazsq) 
#       +log(phi - phiUnifa) + log(phiUnifb - phi)
#       +log(nu - nuUnifa) + log(nuUnifb - nu) 
#       ##Likelihood
#       + lLik
#     )
#     
#     return(logPost)
#   }
#   
#   inits <- c(rep(0, n), starting$beta, rep(0, n), rep(0, n), 
#              log(1), log(1), 
#              logit_bounds(starting$phi, phiUnifa, phiUnifb),
#              logit_bounds(starting$nu, nuUnifa, nuUnifb))
#   
#   chain <- adaptMetropGibbs(ltd = target, starting = inits,
#                             batch = n.batch, 
#                             batch.length = batch.length, 
#                             report = n.report)
#   
#   post_beta <- chain$p.theta.samples[, (n+1):(n+p)]
#   post_z <- chain$p.theta.samples[, (n+p+1):(2*n+p)]
#   post_xi <- chain$p.theta.samples[, 1:n]
#   post_sigmabetasq <- exp(chain$p.theta.samples[, 2*n+p+1])
#   post_sigmazsq <- exp(chain$p.theta.samples[, 2*n+p+2])
#   post_phi <- logitinv_bounds(chain$p.theta.samples[, 2*n+p+3], phiUnifa, phiUnifb)
#   post_nu <- logitinv_bounds(chain$p.theta.samples[, 2*n+p+4], nuUnifa, nuUnifb)
#   
#   return(list(beta = post_beta,
#               z = post_z,
#               xi = post_xi,
#               sigmabetasq = post_sigmabetasq,
#               sigmazsq = post_sigmazsq,
#               phi = post_phi,
#               nu = post_nu))
# }

# library(spBayes)
# library(geoR)

spParams_adaMCMC <- function(z, theta_starting, distmat,
                             n.batch, batch.length, n.report,
                             phiUnifa, phiUnifb, nuUnifa, nuUnifb, nu_z){
  n <- length(as.numeric(z))
  nspParams <- 2
  
  target <- function(theta){
    phi <- logitinv_bounds(theta[1], phiUnifa, phiUnifb)
    nu <- logitinv_bounds(theta[2], nuUnifa, nuUnifb)
    
    Vz <- matern(u = distmat, phi = 1/phi, kappa = nu)
    Lz <- chol(Vz)
    Qz <- sum(backsolve(Lz, z, transpose = TRUE)^2)
    
    logdetz <- 2 * sum(log(diag(Lz)))
    logPost <- -0.5*logdetz - 0.5*(nu_z + n)*log(1 + Qz/nu_z)
    logPost <- logPost + log(phi - phiUnifa) + log(phiUnifb - phi)
    logPost <- logPost + log(nu - nuUnifa) + log(nuUnifb - nu)
    
    return(logPost)
  }
  
  inits <- array(dim = nspParams)
  inits[1] <- logit_bounds(theta_starting[1], phiUnifa, phiUnifb)
  inits[2] <- logit_bounds(theta_starting[2], nuUnifa, nuUnifb)
  
  chain <- adaptMetropGibbs(ltd = target, starting = inits,
                            batch = n.batch, 
                            batch.length = batch.length,
                            verbose = FALSE)
  
  post_phi <- logitinv_bounds(as.numeric(chain$p.theta.samples[, 1]), 
                              phiUnifa, phiUnifb) 
  post_nu <- logitinv_bounds(as.numeric(chain$p.theta.samples[, 2]),
                             nuUnifa, nuUnifb)
  
  return(list(phi = post_phi,
              nu = post_nu))
}

# n <- 100
# coords <- cbind(runif(n), runif(n))
# phi_true <- 3.5
# sigma.sq <- 2
# distmat <- as.matrix(dist(coords))
# R <- sigma.sq * exp(- phi_true * distmat)
# w <- rmvn(1, rep(0,n), R)
# 
# phiUnifa <- 1
# phiUnifb <- 10
# nuUnifa <- 0.1
# nuUnifb <- 2
# nu_z <- 2.1
# 
# out <- spParams_adaMCMC(z = w, theta_starting = c(2, 1), 
#                         distmat = distmat, n.batch = 5, 
#                         batch.length = 10)

# plot(out$phi, type = "l")
# plot(out$nu, type = "l")
# 
# plot(density(out$phi))
# plot(density(out$nu))


spGCM_adaMetropGibbs <- function(y, X, S, family,
                                 N.samp,
                                 spCov,
                                 n_binom = NULL,
                                 marginalPrior = TRUE,
                                 starting,
                                 prior,
                                 Rfastparallel = FALSE){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # prior for spatial process parameters
  phiUnifa <- prior$phi[1]
  phiUnifb <- prior$phi[2]
  nuUnifa <- prior$nu[1]
  nuUnifb <- prior$nu[2]
  
  # prior for variance parameters
  nu_xi <- prior$nu_xi
  nu_beta <- prior$nu_beta
  nu_z <- prior$nu_z
  
  # boundary adjustment parameter
  alpha_epsilon <- prior$alpha_epsilon
  
  if(marginalPrior){
    if(spCov == "exponential"){
      nsp_params <- 1
    }else if(spCov == "matern"){
      nsp_params <- 2
    }else{ stop("Correlation kernel not supported.")}
  }
  
  post_phi <- array(dim = N.samp + 1)
  post_nu <- array(dim = N.samp + 1)
  post_xi <- array(dim = c(n, N.samp + 1))
  post_beta <- array(dim = c(p, N.samp + 1))
  post_z <- array(dim = c(n, N.samp + 1))
  
  post_phi[1] <- starting$phi
  if(spCov == "matern"){
    post_nu[1] <- starting$nu
  }
  post_xi[, 1] <- rep(0, n)
  post_beta[, 1] <- starting$beta
  post_z[, 1] <- rep(0, n)
  
  # Cholesky pre-processing
  distmat <- Rfast::Dist(S)
  if(spCov == "exponential"){
    V_z <- matern(u = distmat, phi = 1/post_phi[1], kappa = 0.5)
  }else if(spCov == "matern"){
    V_z <- matern(u = distmat, phi = 1/post_phi[1], kappa = post_nu[1])
  }else{
    stop("Correlation kernel not supported.")
  }
  L_z <- Rfast::cholesky(V_z, parallel = Rfastparallel)
  XtXplusI <- crossprod(X) / 3 + diag(p)
  XtXplusIchol <- chol(XtXplusI)
  
  n.batch <- 3
  batch.length <- 5
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = N.samp,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = TRUE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  for(s in 1:N.samp){
    
    pb$tick()
    
    # draw posterior of (gamma,q) by Exact Posterior Regression
    
    w_xi <- rnorm(n, mean = 0, sd = nu_xi)
    sigmasq_beta <- 1/rgamma(1, 0.5 * nu_beta, 0.5 * nu_beta)
    sigmasq_z <- 1/rgamma(1, 0.5 * nu_z, 0.5 * nu_z)
    
    if(family == "poisson"){
      w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                         alpha = y[x] + alpha_epsilon, 
                                         kappa = 1, psi = "psi3"))
    }else if(family == "binomial"){
      if(is.null(n_binom)){stop("Supply n_binom.")}
      w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                         alpha = y[x] + alpha_epsilon, 
                                         kappa = n_binom[x] + 2*alpha_epsilon, 
                                         psi = "psi2"))
    }else{ stop("Family currently not supported.") }
    
    w_beta <- rnorm(p) * sqrt(sigmasq_beta)
    w_z <- rnorm(n) * sqrt(sigmasq_z)
    w_z <- as.numeric(crossprod(L_z, w_z))
    
    gamma <- projection(X = X, prechol = XtXplusIchol, 
                        wy = w_y, wxi = w_xi, 
                        wbeta = w_beta, wz = w_z)
    
    post_xi[, s+1] <- gamma[1:n, ]
    post_beta[, s+1] <- gamma[(n+1):(n+p), ]
    post_z[, s+1] <- gamma[(n+p+1):(2*n+p), ]
    
    
    # draw posterior of spatial process parameters
    
    theta_chain <- spParams_adaMCMC(z = post_z[, s+1], 
                                    theta_starting = c(post_phi[s], post_nu[s]), 
                                    distmat = distmat, n.batch = n.batch, 
                                    batch.length = batch.length,
                                    phiUnifa = phiUnifa, 
                                    phiUnifb = phiUnifb, 
                                    nuUnifa = nuUnifa, 
                                    nuUnifb = nuUnifb, 
                                    nu_z = nu_z)
    
    id <- sample(n.batch*batch.length, 1)
    post_phi[s+1] <- theta_chain$phi[id]
    post_nu[s+1] <- theta_chain$nu[id]
  }
  
  return(list(beta = post_beta,
              z = post_z,
              xi = post_xi,
              phi = post_phi,
              nu = post_nu))
}

