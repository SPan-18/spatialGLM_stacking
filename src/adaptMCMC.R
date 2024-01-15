logit_bounds <- function(x, a, b){
  return(log(x - a) - log(b - x))
}

logitinv_bounds <- function(x, a, b){
  return(a + (b - a) * (1 / (1 + exp(- x))))
  # return(b - (b - a) / (1 + exp(x)))
}

spGCM_aMCMC <- function(y, X, S, N.samp, family, 
                        # spCov,
                        n_binom = NULL,
                        marginalPrior = TRUE,
                        starting,
                        prior,
                        tuning,
                        n.batch,
                        batch.length,
                        n.report,
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
  
  H <- rbind(cbind(diag(n), X, diag(n)),
             cbind(diag(n), matrix(0, n, p), matrix(0, n, n)),
             cbind(matrix(0, p, n), diag(p), matrix(0, p, n)),
             cbind(matrix(0, n, n), matrix(0, n, p), diag(n)))
  XtX <- crossprod(X, X)
  HtH <- rbind(cbind(2*diag(n), X, diag(n)),
               cbind(t(X), XtX + diag(p), t(X)),
               cbind(diag(n), X, 2*diag(n)))
  HtHchol <- Rfast::cholesky(HtH, parallel = Rfastparallel)
  PH_sqrt <- backsolve(HtHchol, t(H), transpose = TRUE)
  PH <- crossprod(PH_sqrt)
  QQt <- diag(3*n + p) - PH
  eig_QQt <- eigen(QQt)
  Q <- eig_QQt$vectors[, 1:n]
  
  distmat <- Rfast::Dist(S)
  
  target <- function(theta){
    xi <- theta[1:n]
    beta <- theta[(n+1):(n+p)]
    z <- theta[(n+p+1):(2*n+p)]
    sigmabetasq <- exp(theta[2*n+p+1])
    sigmazsq <- exp(theta[2*n+p+2])
    phi <- logitinv_bounds(theta[2*n+p+3], phiUnifa, phiUnifb)
    nu <- logitinv_bounds(theta[2*n+p+4], nuUnifa, nuUnifb)
    # q <- theta[(2*n+p+5):(3*n+p+4)]
    
    Vz <- matern(u = distmat, phi = 1/phi, kappa = nu)
    Lz <- Rfast::cholesky(Vz, parallel = Rfastparallel)
    Qz <- sum(backsolve(Lz, z, transpose = TRUE)^2)
    Qbeta <- sum(beta^2)
    # Qq <- Q %*% q
    # mu <- - 1.0 * Qq[1:n]
    logLambda <- X %*% beta + z + xi #- mu
    if(family == "poisson"){
      Qxi <- sum(xi^2)
      lLik <- sum(dpois(y, exp(logLambda), log = TRUE))
    }
    
    logPost <- (
      ##Priors
      -(0.5 * nu_beta + 1) * log(sigmabetasq) - (0.5 * nu_beta)/sigmabetasq
      -(0.5 * nu_z + 1) * log(sigmazsq) - (0.5 * nu_z)/sigmazsq
      - 0.5 * n * log(sigmazsq) - 0.5 * sum(log(diag(Lz))) - 0.5 * Qz/sigmazsq
      - 0.5 * p * log(sigmabetasq) - 0.5 * Qbeta/sigmabetasq
      # + sum(alpha_epsilon * logLambda) 
      - 0.5 * Qxi / (nu_xi)^2
      ##Jacobians
      +log(sigmabetasq) + log(sigmazsq) 
      +log(phi - phiUnifa) + log(phiUnifb - phi)
      +log(nu - nuUnifa) + log(nuUnifb - nu) 
      ##Likelihood
      + lLik
    )
    
    return(logPost)
  }
  
  inits <- c(rep(0, n), starting$beta, rep(0, n), rep(0, n), 
             log(1), log(1), 
             logit_bounds(starting$phi, phiUnifa, phiUnifb),
             logit_bounds(starting$nu, nuUnifa, nuUnifb))
  
  chain <- adaptMetropGibbs(ltd = target, starting = inits,
                            batch = n.batch, 
                            batch.length = batch.length, 
                            report = n.report)
  
  post_beta <- chain$p.theta.samples[, (n+1):(n+p)]
  post_z <- chain$p.theta.samples[, (n+p+1):(2*n+p)]
  post_xi <- chain$p.theta.samples[, 1:n]
  post_sigmabetasq <- exp(chain$p.theta.samples[, 2*n+p+1])
  post_sigmazsq <- exp(chain$p.theta.samples[, 2*n+p+2])
  post_phi <- logitinv_bounds(chain$p.theta.samples[, 2*n+p+3], phiUnifa, phiUnifb)
  post_nu <- logitinv_bounds(chain$p.theta.samples[, 2*n+p+4], nuUnifa, nuUnifb)
  
  return(list(beta = post_beta,
              z = post_z,
              xi = post_xi,
              sigmabetasq = post_sigmabetasq,
              sigmazsq = post_sigmazsq,
              phi = post_phi,
              nu = post_nu))
}

