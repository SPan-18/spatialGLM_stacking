n <- 100
p <- 2
Nt <- 3
time <- rep(1:Nt, each = n)
X <- cbind(rep(1, n*Nt), sapply(1:(p-1), function(x) rnorm(n*Nt)))
beta0 <- c(4,5,6,-0.5,0.5,1.5)
bigX <- mkdynX(X, time)
Xbeta <- bigX %*% beta0
y <- Xbeta + rnorm(n*Nt)

XtX <- crossprod(bigX)
Xty <- crossprod(bigX, y)
beta_hat <- mysolve(XtX, Xty)
print(beta_hat)
