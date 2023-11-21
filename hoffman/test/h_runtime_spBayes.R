rm(list = ls())

library(spBayes)

dat <- read.csv("../data/h_sim_count5000.csv")

nseq <- c(c(1:10)*100, 2:5*1000)
n_run <- length(nseq)
runtime <- array(dim = c(n_run, 2))
runtime[, 1] <- nseq
colnames(runtime) <- c("n", "time")

for(i in 1:n_run){
  simdat <- dat[1:nseq[i], ]
  y <- as.numeric(simdat$y)
  X <- as.matrix(simdat[, grep("x", names(simdat))])
  X <- as.matrix(X[, -1])
  S <- as.matrix(simdat[, c("s1", "s2")])
  
  ##Collect samples
  beta.starting <- coefficients(glm(y~X, family="poisson"))
  beta.tuning <- t(chol(vcov(glm(y~X, family="poisson"))))
  n.batch <- 500
  batch.length <- 50
  n.samples <- n.batch*batch.length
  
  cat("Running n =", nseq[i], "...\t")
  t1 <- Sys.time()
  m.1 <- spGLM(y~X, family="poisson", coords=S,
               starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
               tuning=list("beta"=c(0.1, 0.1), "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
               priors=list("beta.Flat", "phi.Unif"=c(0.03, 10), "sigma.sq.IG"=c(2, 1)),
               amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
               cov.model="exponential", verbose=FALSE, n.report=10)
  t2 <- Sys.time()
  rt <- difftime(t2, t1, units = "secs")
  runtime[i, 2] <- rt
  cat("Took", rt, units(rt), "\n")
}

write.csv(runtime, "spBayes_runtime.csv", row.names = FALSE)
