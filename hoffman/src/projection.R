# fast computation of inv(H'H)H'y

projection <- function(X, prechol, wy, wxi, wbeta, wz){
  p <- dim(X)[2]
  w1_tilde <- 0.5 * wy + wxi - 0.5 * wz
  w2_tilde <- 0.5 * crossprod(X, wy) + wbeta - 0.5 * crossprod(X, wz)
  w3_tilde <- (1/3) * crossprod(X, w1_tilde) - w2_tilde
  if(is.null(prechol)){
    temp_pp <- crossprod(X) / 3
    diag(temp_pp) <- diag(temp_pp) + 1
    w3_tilde <- mysolve(temp_pp, w3_tilde)
  }else{
    w3_tilde <- backsolve(prechol, w3_tilde, transpose = TRUE)
    w3_tilde <- backsolve(prechol, w3_tilde)
  }
  temp <- X %*% w3_tilde
  w4_tilde <- (2/3) * w1_tilde + (1/3) * temp
  temp2 <- - 0.5 * (w4_tilde - temp) + 0.5 * (wy + wz)
  return(rbind(w4_tilde, - w3_tilde, temp2))
}

# TEST CORRECTNESS:
# n = 3
# p = 2
# X = cbind(rep(1,n), 1:n)
# wy = rep(1, n)
# wxi = rep(1, n)
# wbeta = rep(1, p)
# wz = rep(1, n)
# 
# HtH <- rbind(cbind(2*diag(n), X, diag(n)),
#              cbind(t(X), crossprod(X) + diag(p), t(X)),
#              cbind(diag(n), X, 2*diag(n)))
# Htw <- c(wy + wxi,
#          as.numeric(crossprod(X, wy)) + wbeta,
#          wy + wz)
# 
# XtXplusI <- crossprod(X) / 3 + diag(p)
# XtXplusIchol <- chol(XtXplusI)
# mysolve(HtH, Htw)
# projection(X, prechol = NULL, wy, wxi, wbeta, wz)
# projection(X, prechol = XtXplusIchol, wy, wxi, wbeta, wz)
