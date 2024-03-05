ci_beta <- function(mat){
  out <- array(dim = c(ncol(mat), 5))
  for(i in 1:ncol(mat)){
    out[i, ] = c(mean(mat[, i]), sd(mat[, i]), 
                 quantile(mat[, i], probs = c(0.025, 0.5, 0.975)))
  }
  colnames(out) <- c("Mean", "SD", "2.5%", "50%", "97.5%")
  rownames(out) <- paste("beta", 1:ncol(mat) - 1, sep = '')
  return(out)
}

logit.inv <- function(x) return(exp(x)/(1 + exp(x)))

kbar_binary <- function(alpha, kappa){
  return(digamma(alpha) - digamma(kappa - alpha) * beta(alpha, kappa - alpha))
}

# solve Ax = b, A symmetric p.d.
mysolve <- function(A, b){
  a <- chol(A)
  x <- backsolve(a, b, transpose = TRUE)
  x <- backsolve(a, x)
  return(x)
}

# divide a set into K equal parts
id_partition <- function(N, K, random = TRUE){
  if(random){
    ids <- 1:N
    idlist <- split(ids, sort(ids %% K))
    permt <- sample(ids)
    for(i in 1:K){
      idlist[[i]] <- permt[idlist[[i]]]
    }
    names(idlist) = NULL
  }else{
    ids <- 1:N
    idlist <- split(ids, sort(ids %% K))
    names(idlist) = NULL
  }
  return(idlist)
}

# G_phi <- c(2, 3)
# G_alphaepsilon <- c(0.5, 0.75)
# model_list <- create_model_list(G_phi, G_alphaepsilon)

create_model_list <- function(G_decay, G_smoothness, G_timedecay = NULL,
                              G_rho = NULL,
                              G_epsilon, G_nuxi = 0,
                              G_nubeta = 1, G_nuz = 1){
  
  if(is.null(G_timedecay)){
    models <- expand.grid(G_decay, G_smoothness, G_epsilon, G_nuxi, G_nubeta, G_nuz)
    names(models) <- c("phi", "smooth", "alpha_epsilon", "nu_xi", "nu_beta", "nu_z")
    model_list <- lapply(1:nrow(models), function(x){
      return(list(phi = models[x, "phi"],
                  nu_matern = models[x, "smooth"],
                  alpha_epsilon = models[x, "alpha_epsilon"],
                  nu_xi = models[x, "nu_xi"],
                  nu_beta = models[x, "nu_beta"],
                  nu_z = models[x, "nu_z"]))})
  }else{
    models <- expand.grid(G_decay, G_smoothness, G_timedecay, G_rho,
                          G_epsilon, G_nuxi, G_nubeta, G_nuz)
    names(models) <- c("phi", "smooth", "phi_t", "rho", "alpha_epsilon", 
                       "nu_xi", "nu_beta", "nu_z")
    model_list <- lapply(1:nrow(models), function(x){
      return(list(phi = models[x, "phi"],
                  nu_matern = models[x, "smooth"],
                  phi_t = models[x, "phi_t"],
                  rho = models[x, "rho"],
                  alpha_epsilon = models[x, "alpha_epsilon"],
                  nu_xi = models[x, "nu_xi"],
                  nu_beta = models[x, "nu_beta"],
                  nu_z = models[x, "nu_z"]))})
  }
  
  return(model_list)
}

mkdynX <- function(X, time){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  Nt <- max(time)
  
  X_big <- matrix(0, n, p * Nt)
  
  for(i in 1:Nt){
    ids <- which(time == i)
    for(j in 1:p){
      X_big[ids, (j-1)*Nt + i] <- X[ids, j]
    }
  }
  
  return(X_big)
}

mkdynX2 <- function(X, time){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  Nt <- max(time)
  
  X_big <- matrix(0, n, p * Nt)
  
  for(i in 1:Nt){
    ids <- which(time == i)
    X_big[ids, ((i-1)*p+1):(i*p)] <- X[ids,]
  }
  
  return(X_big)
}

# nX <- 6
# X1 <- cbind(rep(1, nX), rnorm(nX), rnorm(nX))
# time_1 <- rep(1:3, each = 2)
# mkdynX2(X1, time_1)


mkbetacov <- function(size, rho, prior = "AR1"){
  if(prior == "AR1"){
    out <- as.matrix(dist(1:size))
    out <- rho^out #/ (1-rho^2)
  }else stop("Invalid prior for covariance structure of beta.")
  return(out)
}

# mkbetacov(4, 0.5)

makeG <- function(X){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  G <- matrix(0, n, n * p)
  
  for(i in 1:p){
    ids <- ((i-1)*n + 1):(i*n)
    G[, ids] <- diag(X[, i])
  }
  
  return(G)
}

# nX <- 5
# X1 <- cbind(rep(1, nX), rnorm(nX), rnorm(nX))
# makeG(X1)

# XtG <- function(X, tr = F){
#   n <- dim(X)[1]
#   p <- dim(X)[2]
#   out <- matrix(0, p, n*p)
#   for(i in 1:p){
#     ids <- ((i-1)*n+1):(i*n)
#     for(j in 1:n){
#       out[, ids[j]] <- X[j, ] * X[j, i]
#     }
#   }
#   if(tr){
#     return(t(out))
#   }else{
#     return(out)
#   }
# }

# nX <- 5
# X1 <- cbind(rep(1, nX), rnorm(nX), rnorm(nX))
# all(crossprod(X1, makeG(X1)) == XtG(X1))

GtGplusI <- function(X){
  G <- makeG(X)
  GtG <- crossprod(G)
  diag(GtG) <- diag(GtG) + 1
  return(GtG)
}
