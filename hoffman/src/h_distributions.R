# functions related to the Diaconis-Ylvisaker distribution
# see Diaconis, Ylvisaker (1979); DOI: 10.1214/aos/1176344611
# notations follow Bradley et. al. (2020)

# The following is a conjugate prior for a data model from natural
# exponential family; alpha, kappa: DY parameters; K(alpha): normalizing constant
# f(y|alpha, kappa) = K(alpha, kappa) exp{alpha*Y - kappa*psi(y)}, ...
# dDY gives the value of the density of DY distribution
# rDY generates random samples from the DY distribution

# Choices of unit log-partition function of data model considered:
# 1. psi1(y) = log(-1/y)      (data model: gamma)
# 2. psi2(y) = log(1+exp(y))  (data model: binomial)
# 3. psi3(y) = exp(y)         (data model: poisson)
# 4. psi4(y) = y^2            (data model: normal)

# Unit log-partition functions
psi1 <- function(y){
  if(all(y) < 0) return(log(- (1 / y)))
  else stop("Invalid prior support: psi1")
}

psi2 <- function(y) log(1 + exp(y))

psi3 <- function(y) exp(y)

psi4 <- function(y) y^2

# random sample from DY distribution with parameters alpha and kappa 
# with unit log-partition function psi
rDY <- function(n.samples, alpha, kappa, psi){
  if(psi == "psi1"){
    if(alpha > 0 && kappa > 0){
      w <- rgamma(n = n.samples, shape = kappa + 1, rate = alpha)
      return(- w)
    }else stop("Invalid DY parameters input (psi1).")
  }else if(psi == "psi2"){
    if(kappa > alpha && alpha > 0){
      w <- rbeta(n = n.samples, shape1 = alpha, shape2 = kappa - alpha)
      return(log(w / (1 - w)))
    }else{
      stop("Invalid DY parameters input (psi2).")
    } 
  }else if(psi == "psi3"){
    if(alpha > 0 && kappa > 0){
      w <- rgamma(n = n.samples, shape = alpha, rate = kappa)
      return(log(w))
    }else stop("Invalid DY parameters input (psi3).")
  }else if(psi == "psi4"){
    if(kappa > 0){
      w <- rnorm(n = n.samples, mean = alpha/(2*kappa), 
                 sd = sqrt(1/(2*kappa)))
      return(w)
    }else stop("Invalid DY parameters input (psi4).")
  }else stop("Undefined unit log-partition function (psi).")
}

# # test rDY; high alpha and kappa will result in asymptotic normality
# rdat <- rDY(n.samples = 1000, alpha = 100, kappa = 150, psi = "psi1")
# hist(rdat)

# normalizing constant K(alpha, kappa)

normconstDY <- function(alpha, kappa, psi){
  if(psi == "psi1"){
    if(alpha > 0 && kappa > 0){
      k <- alpha^(kappa + 1) / gamma(kappa + 1)
      return(k)
    }else stop("Invalid DY parameters input (psi1).")
  }else if(psi == "psi2"){
    if(kappa > alpha && alpha > 0){
      k <- gamma(kappa) / (gamma(alpha) * gamma(kappa - alpha))
      return(k)
    }else stop("Invalid DY parameters input (psi2).")
  }else if(psi == "psi3"){
    if(alpha > 0 && kappa > 0){
      k <- kappa^alpha / gamma(alpha)
      return(k)
    }else stop("Invalid DY parameters input (psi3).")
  }else if(psi == "psi4"){
    if(kappa > 0){
      k <- sqrt(kappa / pi) * exp(- alpha^2 / (4 * kappa))
      return(k)
    }else stop("Invalid DY parameters input (psi4).")
  }else stop("Undefined unit log-partition function (psi).")
}

# density of a DY random variable

dDY <- function(x, alpha, kappa, psi){
  K <- normconstDY(alpha, kappa, psi)
  psifn <- get(psi)
  out <- (alpha * x) - (kappa * psifn(x))
  return(K * exp(out))
}

# normdat <- rnorm(10)
# dn <- dnorm(normdat)
# dyn <- dDY(normdat, 0, 0.5, "psi4")

# functions related to the CM (Conjugate Multivariate) distribution
# multivariate extension of DY random variable
# see Bradley, Holan, Wikle (2020); DOI: 10.1080/01621459.2019.1677471
# y ~ CM if y = m + Lw, where w_i ~ DY

rCM <- function(n.samples, mu, V, alpha, kappa, psi){
  d <- length(alpha)
  temp <- array(dim = c(d, n.samples))
  for(i in 1:d){
    temp[i, ] <- rDY(n.samples, alpha[i], kappa[i], psi)
  }
  temp <- V %*% temp
  temp <- apply(temp, 2, function(x) x + mu)
  return(t(temp))
}

# rCM(10, c(0, 1), V = matrix(1, 2, 2), alpha = c(0, 0),
#     kappa = c(1, 2), psi = "psi4")

# normalizing constant det(V.inv) * Prod_i K(alpha_i, kappa_i)
normconstCM <- function(alpha, kappa, V, psivec, log = FALSE){
  d <- length(alpha)
  k <- array(dim = d)
  for(i in 1:d){
    k[i] <- normconstDY(alpha[i], kappa[i], psivec[i])
  }
  if(is.lower.tri(V, diag = TRUE)){
    logdet <- log(prod(diag(V)))
  }else{
    logdet <- log(det(V))
  }
  logconst <- - logdet + sum(log(k))
  if(log) return(logconst)
  else return(exp(logconst))
}

# calculate density of CM, needs V to be lower triangular
dCM <- function(x, mu, V, alpha, kappa, psivec){
  K <- normconstCM(alpha, kappa, V, psivec)
  if(is.vector(x)){
    y <- forwardsolve(V, x - mu)
    psi.y <- array(dim = length(psivec))
    for(i in 1:length(psivec)){
      psifn <- get(psivec[i])
      psi.y[i] <- psifn(y[i])
    }
    temp <- crossprod(alpha, y) - crossprod(kappa, psi.y)
  }else if(is.matrix(x)){
    n <- nrow(x)
    y.temp <- t(x) - tcrossprod(mu, rep(1, n))
    y <- t(forwardsolve(V, y.temp))
    temp <- array(dim = n)
    for(i in 1:n){
      psi.y <- array(dim = length(psivec))
      for(j in 1:length(psivec)){
        psifn <- get(psivec[j])
        psi.y[j] <- psifn(y[i, j])
      }
      temp[i] <- crossprod(alpha, y[i, ]) - crossprod(kappa, psi.y)
    }
  }
  return(as.numeric(K * exp(temp)))
}

# V.trial <- matrix(c(1, 0, 0, 2, 1, 0, 2, 3, 1), nrow = 3, byrow = TRUE)
# normmat <- mvrnorm(n = 1, mu = c(1, 2, 0), Sigma = tcrossprod(V.trial))
# psi.vec <- rep("psi4", 3)
# library(emdbook)
# dmvnorm(normmat, c(0,0,0), diag(c(1,4,1)))
# dCM(normmat, c(0,0,0), diag(c(1, 2, 1)), rep(0, 3), rep(0.5, 3), psi.vec)
# psi.val <- "psi2"
# psi.vec <- rep(psi.val, 3)
# cmdat <- rCM(10, c(1,1,1), diag(3), c(1,2,5), rep(20,3), psi.val)
# dCM(cmdat, rep(0,3), diag(3), rep(5,3), rep(15,3), psi.vec)


# Poisson density with b ≠ 1
dpois_b <- function(x, lambda, b, log = FALSE){
  logadj <- (x * log(b)) - ((b - 1) * lambda)
  if(log) return(dpois(x, lambda, log = TRUE) * logadj)
  else return(dpois(x, lambda, log = FALSE) * exp(logadj))
}

# xseq <- 1:20
# y1 <- dpois_b(xseq, 5, 1)
# y2 <- dpois_b(xseq, 5, 0.5)
# y3 <- dpois_b(xseq, 5, 1.5)
# 
# plot(xseq, y1, col = "red", type = "o", ylim = c(0, 0.3), cex = 0.8,
#      xlab = bquote(y), ylab = bquote(p~"("~y~"|"~b~","~lambda~")"))
# lines(xseq, y2, col = "blue", type = "o", cex = 0.8,)
# lines(xseq, y3, col = "darkgreen", type = "o", cex = 0.8,)
# legend(15.5, 0.33, legend = c(bquote(b == 0.5), 
#                               bquote(b == "1.0"), 
#                               bquote(b == 1.5)),
#        col = c("blue", "red", "darkgreen"), lty = c(1, 1, 1),
#        bty = "n", cex = 0.75, y.intersp = 0.5)

