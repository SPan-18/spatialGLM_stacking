# Fast O(nr^3+p^3) evaluation of projection step

trsolve_prechol <- function(cholA, b){
  x <- backsolve(cholA, b, transpose = TRUE)
  x <- backsolve(cholA, x)
  return(x)
}

projection2 <- function(X, X_tilde, chol_p, chol_nr, 
                        v_eta, v_xi, v_beta, v_z){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  r <- dim(X_tilde)[2]
  
  v2 <- sapply(1:n, function(x){
    (v_eta[x] * X_tilde[x, ]) + v_z[((x-1)*r+1):(x*r)] })
  
  v3 <- sapply(1:n, function(x){
    trsolve_prechol(chol_nr[[x]], v2[, x])
  })
  
  v4 <- sapply(1:n, function(x){ crossprod(X_tilde[x, ], v3[, x]) })
  
  v51 <- v_eta - v4
  v52 <- crossprod(X, v51) + v_beta
  v51 <- v51 + v_xi
  
  diag1 <- apply(X_tilde, 1, function(x){ 1/(2 + sum(x^2)) })
  v6 <- crossprod(X, diag1 * v51) - v52
  
  v7 <- trsolve_prechol(chol_p, v6)
  
  diag2 <- apply(X_tilde, 1, function(x){ (1 + sum(x^2))/(2 + sum(x^2)) })
  v8 <- diag2 * v51 + diag1 * (X %*% v7)
  
  v9 <- v8 - (X %*% v7)
  
  v10 <- sapply(1:n, function(x){ 
    trsolve_prechol(chol_nr[[x]], v9[x]*X_tilde[x, ]) })
  v10 <- - v10 + v3
  
  return(c(as.numeric(v8), - as.numeric(v7),
           as.numeric(c(v10))))
}

# test projection2
# set.seed(1729)
# n <- 1000
# p <- 3
# r <- 3
# X1 <- cbind(rep(1, n), sapply(1:(p-1), function (x) rnorm(n)))
# X_tilde_1 <- X1
# v_eta_1 <- rnorm(n)
# v_xi_1 <- rnorm(n)
# v_beta_1 <- rnorm(p)
# v_z_1 <- rnorm(n*r)
# 
# diag1 <- apply(X_tilde_1, 1, function(x){ 1/(2 + sum(x^2)) })
# S_D_star <- crossprod(X1, diag1 * X1)
# diag(S_D_star) <- diag(S_D_star) + 1
# chol_p_1 <- chol(S_D_star)
# chol_nr_1 <- lapply(1:n, function(x){chol(tcrossprod(X_tilde_1[x, ]) + diag(r))})
# 
# res <- projection2(X = X1, X_tilde = X_tilde_1,
#                    chol_p = chol_p_1, chol_nr = chol_nr_1,
#                    v_eta = v_eta_1, v_xi = v_xi_1,
#                    v_beta = v_beta_1, v_z = v_z_1)
# 
# X_tilde_big <- matrix(0, n, n*r)
# for(i in 1:n){
#   X_tilde_big[i, ((i-1)*r + 1):(i*r)] = X_tilde_1[i, ]
# }
# res2 <- X_tilde_big %*% solve(crossprod(X_tilde_big) + diag(n*r)) %*% (crossprod(X_tilde_big, v_eta_1) + v_z_1)
# H_big <- cbind(diag(n), X1, X_tilde_big)
# res2 <- solve(crossprod(H_big) + diag(n+p+n*r)) %*%
#   (crossprod(H_big, v_eta_1) + c(v_xi_1, v_beta_1, v_z_1))
# summary(as.numeric(c(res2)) - res)

# P_tilde <- X_tilde_big %*% solve(crossprod(X_tilde_big) + diag(n*r)) %*% t(X_tilde_big)
# # (diag(P_tilde) - apply(X_tilde_1, 1, function(x){ sum(x^2)/(1 + sum(x^2)) }))
# S_d_star2 <- t(X1) %*% (diag(n) - P_tilde) %*% X1 + diag(p) -
#   t(X1) %*% (diag(n) - P_tilde) %*% solve(diag(2, n) - P_tilde) %*%
#   (diag(n) - P_tilde) %*% X1


