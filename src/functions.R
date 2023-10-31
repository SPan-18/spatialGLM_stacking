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
id_partition <- function(N, K){
  n <- floor(N / K)
  boxid <- c(rep(1:K, n), rep(1, N %% K))
  boxes <- sample(boxid)
  idlist <- vector(mode = "list", length = K)
  for(i in 1:K){
    idlist[[i]] <- which(boxes == i)
  }
  return(idlist)
}

create_model_list <- function(G_phi, G_alphaepsilon, G_nuxi = 0,
                              G_nubeta = 0.01, G_nuz = 0.01){
  models <- expand.grid(G_phi, G_alphaepsilon, G_nuxi, G_nubeta, G_nuz)
  names(models) <- c("phi", "alpha_epsilon", "nu_xi", "nu_beta", "nu_z")
  model_list <- lapply(1:nrow(models), function(x){
    return(list(phi = models[x, "phi"],
                alpha_epsilon = models[x, "alpha_epsilon"],
                nu_xi = models[x, "nu_xi"],
                nu_beta = models[x, "nu_beta"],
                nu_z = models[x, "nu_z"]))})
  return(model_list)
}
