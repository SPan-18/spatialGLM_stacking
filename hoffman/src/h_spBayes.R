rm(list = ls())

library(MBA)
library(coda)
library(spBayes)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p) %*% D + rep(mu,rep(n,p)))
}

###########################
##Spatial poisson
###########################
##Generate count data
set.seed(1)
n <- 100
coords <- cbind(runif(n,1,100),runif(n,1,100))
phi <- 3/50
sigma.sq <- 2
R <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- rmvn(1, rep(0,n), R)
x <- as.matrix(rep(1,n))
beta <- 0.1
y <- rpois(n, exp(x%*%beta+w))

##Collect samples
beta.starting <- coefficients(glm(y~x-1, family="poisson"))
beta.tuning <- t(chol(vcov(glm(y~x-1, family="poisson"))))
n.batch <- 500
batch.length <- 50
n.samples <- n.batch*batch.length
##Note tuning list is now optional
t1 <- Sys.time()
m.1 <- spGLM(y~1, family="poisson", coords=coords,
             starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=0.1, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Flat", "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)
t2 <- Sys.time()
print(t2-t1)

##Just for fun check out the progression of the acceptance
##as it moves to 43% (same can be seen for the random spatial effects).
plot(mcmc(t(m.1$acceptance)), density=FALSE, smooth=FALSE)

##Now parameter summaries, etc.
burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples
m.1$p.beta.theta.samples[,"phi"] <- 3/m.1$p.beta.theta.samples[,"phi"]
plot(m.1$p.beta.theta.samples)
print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))
beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]
y.hat <- apply(exp(x%*%beta.hat+w.hat), 2, function(x){rpois(n, x)})
y.hat.mu <- apply(y.hat, 1, mean)
##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords,y),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed counts")
contour(surf, add=TRUE)
text(coords, labels=y, cex=1)
surf <- mba.surf(cbind(coords,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Fitted counts")
contour(surf, add=TRUE)
text(coords, labels=round(y.hat.mu,0), cex=1)
