# find Cholesky with block row-column removed

cholesky_CV <- function(L_full, ids){
  n <- dim(L_full)[1]
  n_k <- length(ids)
  if(ids[n_k] == n){
    chol_k <- L_full[1:(ids[1] - 1), 1:(ids[1] - 1)]
  }else if(ids[1] == 1){
    L33 <- L_full[(ids[n_k]+1):n, (ids[n_k]+1):n]
    L23 <- L_full[ids, (ids[n_k] + 1):n]
    # chol_k <- chol_rankupdate(L33, L23)
    chol_k <- Rfast::cholesky(crossprod(L33) + crossprod(L23))
  }else{
    chol_k <- matrix(0, n - n_k, n - n_k)
    chol_k[1:(ids[1] - 1), 1:(ids[1] - 1)] = L_full[1:(ids[1] - 1), 1:(ids[1] - 1)]
    chol_k[1:(ids[1] - 1), (ids[1]:(n - n_k))] = 
      L_full[1:(ids[1] - 1), (ids[n_k]+1):n]
    L33 <- L_full[(ids[n_k]+1):n, (ids[n_k]+1):n]
    L23 <- L_full[ids, (ids[n_k] + 1):n]
    # chol_k[(ids[1]:(n - n_k)), (ids[1]:(n - n_k))] = chol_rankupdate(L33, L23)
    chol_k[(ids[1]:(n - n_k)), (ids[1]:(n - n_k))] <- Rfast::cholesky(crossprod(L33) + crossprod(L23))
  }
  return(chol_k)
}

# finds Cholesky of A + M'M with rank updates
### SLOW! Can we bypass `for` loop?
chol_rankupdate <- function(cholA, M){
  for(i in 1:nrow(M)){
    cholA <- cholupdate(cholA, M[i, ])
  }
  return(cholA)
}

# library(fastmatrix)

# a <- matrix(c(1,1,1,1,2,3,1,3,6), ncol = 3)
# r <- chol(a)
# x1 <- c(0,0,1)
# x2 <- c(0, 1, 0)
# b <- a + outer(x1,x1) + outer(x2, x2)
# r1 <- cholupdate(r, x1)
# r2 <- cholupdate(r1, x2)
# r2
# all(r2 == chol(b))
# all(chol(b) == chol_rankupdate(r, rbind(x1, x2)))


# n <- 1000
# A <- matrix(rnorm(n^2), n, n)
# A <- crossprod(A)
# M <- matrix(rnorm(0.5*n^2), n/2, n)
# cholA <- Rfast::cholesky(A)
# t1 <- Sys.time()
# res1 <- Rfast::cholesky(A + crossprod(M))
# t2 <- Sys.time()
# t2 - t1
# 
# t3 <- Sys.time()
# res2 <- chol_rankupdate(cholA = cholA, M = M)
# t4 <- Sys.time()
# t4 - t3

# TEST cholesky_CV
# set.seed(1729)
# M1 <- matrix(rnorm(49), 7, 7)
# M <- crossprod(M1)
# Mchol <- chol(M)
# ids <- c(1,2,3,4)
# Msub <- M[-ids, -ids]
# Msubchol <- chol(Msub)
# Msubchol2 <- cholesky_CV(Mchol, ids = ids)
# all(abs(Msubchol-Msubchol2) < 1E-15)

# rm(list = ls())
# 
# library(Rfast)

# n <- 7000
# Grid <- data.frame(x = runif(n, 0, 1),
#                    y = runif(n, 0, 1))
# 
# distance <- as.matrix(dist(Grid))
# phi <- 3
# Rphi <- exp(- phi * distance)
# 
# # t11 <- Sys.time()
# # chol_base <- chol(Rphi)
# # t12 <- Sys.time()
# # cat("Base:\n")
# # print(t12 - t11)
# 
# t21 <- Sys.time()
# chol_Rfast <- Rfast::cholesky(Rphi, parallel = TRUE)
# t22 <- Sys.time()
# cat("Rfast (parallel = TRUE):\n")
# print(t22 - t21)
# 
# t31 <- Sys.time()
# chol_Rfast2 <- Rfast::cholesky(Rphi, parallel = FALSE)
# t32 <- Sys.time()
# cat("Rfast (parallel = FALSE):\n")
# print(t32 - t31)
