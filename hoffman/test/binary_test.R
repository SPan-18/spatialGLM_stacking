slow_posterior_GCM <- function(y, X, S, N.samp, 
                               family = "binomial",
                               fixed_prior = "gaussian", 
                               spatial_prior = "gaussian", 
                               mod_params){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  phi <- mod_params$phi
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  nu_z <- mod_params$nu_z
  alpha_epsilon <- mod_params$alpha_epsilon
  # kappa_epsilon <- mod_params$kappa_epsilon = 0
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = N.samp,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  
  distmat <- as.matrix(dist(S))
  XtX <- crossprod(X, X)
  HtH <- rbind(cbind(2*diag(n), X, diag(n)),
               cbind(t(X), XtX + diag(p), t(X)),
               cbind(diag(n), X, 2*diag(n)))
  chol11 <- Sys.time()
  HtHchol <- chol(HtH)
  chol12 <- Sys.time()
  chol21 <- Sys.time()
  L_z <- chol(exp(- phi * distmat))
  chol22 <- Sys.time()
  cat("Cholesky runtimes: \n", 
      "HtH:", chol12-chol11, units(chol12-chol11), "\n",
      "R(\U03C6):", chol22-chol21, units(chol22-chol21), "\n")
  beta_samps <- array(dim = c(N.samp, p))
  z_samps <- array(dim = c(N.samp, n))
  
  for(i in 1:N.samp){
    pb$tick()
    if(nu_xi == 0){
      sigmasq_xi <- 1
    } else{
      sigmasq_xi <- 1/rgamma(1, nu_xi/2, nu_xi/2)
    } 
    sigmasq_beta <- 1/rgamma(1, nu_beta/2, nu_beta/2)
    sigmasq_z <- 1/rgamma(1, nu_z/2, nu_z/2)
    
    w_y <- sapply(1:n, function(x) rDY(n.samples = 1, 
                                       alpha = y[x] + alpha_epsilon, 
                                       kappa = 1 + 2*alpha_epsilon, 
                                       psi = "psi2"))
    
    w_xi <- rnorm(n = n, mean = 0, sd = sqrt(sigmasq_xi))
    w_beta <- rnorm(n = p, mean = 0, sd = sqrt(sigmasq_beta))
    w_z <- rnorm(n = n, mean = 0, sd = sqrt(sigmasq_z))
    w_z <- as.numeric(crossprod(L_z, w_z))
    
    Htw <- c(w_y + w_xi, 
             as.numeric(crossprod(X, w_y)) + w_beta,
             w_y + w_z)
    temp <- forwardsolve(t(HtHchol), Htw)
    gamma <- backsolve(HtHchol, temp)
    beta_samps[i, ] <- gamma[(n+1):(n+p)]
    z_samps[i, ] <- gamma[(n+p+1):(2*n+p)]
  }
  return(list(beta = beta_samps,
              z = z_samps))
}

binomial_GCM_nonsp <- function(y, X, N.samp, mod_params){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  nu_xi <- mod_params$nu_xi
  nu_beta <- mod_params$nu_beta
  alpha_epsilon <- mod_params$alpha_epsilon
  XtX <- crossprod(X, X)
  HtH <- rbind(cbind(2*diag(n), X),
               cbind(t(X), XtX + diag(p)))
  HtHchol <- chol(HtH)
  beta_samps <- array(dim = c(N.samp, p))
  
  for(i in 1:N.samp){
    if(nu_xi == 0){
      sigmasq_xi <- 1
    } else{
      sigmasq_xi <- 1/rgamma(1, nu_xi/2, nu_xi/2)
    } 
    sigmasq_beta <- 1/rgamma(1, nu_beta/2, nu_beta/2)
    w_y <- array(dim = n)
    for(j in 1:n){
      w <- rbeta(1, y[j] + alpha_epsilon, 20 - y[j] + alpha_epsilon)
      w_y[j] <- log(w/(1-w))
    }
    w_xi <- rnorm(n = n, mean = 0, sd = sqrt(sigmasq_xi))
    w_beta <- rnorm(n = p, mean = 0, sd = sqrt(sigmasq_beta))
    Htw <- c(w_y + w_xi, 
             as.numeric(crossprod(X, w_y)) + w_beta)
    temp <- forwardsolve(t(HtHchol), Htw)
    gamma <- backsolve(HtHchol, temp)
    # if(any(gamma[(n+1):(n+p)] == "NaN")) browser()
    beta_samps[i, ] <- gamma[(n+1):(n+p)]
  }
  return(beta_samps)
}

source("../src/distributions.R")
source("../src/functions.R")
dat <- read.csv("../data/sim_nspbinom1000.csv")
# dat <- dat[1:100, ]
# glm_out <- glm(y ~ x1, family = binomial(link = "logit"), data = dat)
# summary(glm_out)
y <- as.numeric(dat$y)
X <- as.matrix(dat[, grep("x", names(dat))])
# S <- as.matrix(dat[, c("s1", "s2")])
mod_params1 <- list(nu_xi = 0, nu_beta = 2,
                    nu_z = 2, alpha_epsilon = 0.5)

out <- binomial_GCM_nonsp(y = y, X = X, N.samp = 100, 
                          mod_params = mod_params1)
beta_samps <- out
print(ci_beta(beta_samps))
