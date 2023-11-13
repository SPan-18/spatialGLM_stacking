library(spBayes)

dat <- read.csv("../data/h_sim_count1000.csv")

n_run <- 2
runtime <- array(dim = c(n_run, 2))
nseq <- c(1:n_run)*100
runtime[, 1] <- nseq
colnames(runtime) <- c("n", "time")

for(i in 1:n_run){
  simdat <- dat[nseq[i], ]
  y <- as.numeric(simdat$y)
  X <- as.matrix(simdat[, grep("x", names(simdat))])
  X <- X[, -1]
  S <- as.matrix(simdat[, c("s1", "s2")])
  
  ##Collect samples
  beta.starting <- coefficients(glm(y~X-1, family="poisson"))
  beta.tuning <- t(chol(vcov(glm(y~X-1, family="poisson"))))
  n.batch <- 500
  batch.length <- 50
  n.samples <- n.batch*batch.length
  t1 <- Sys.time()
  m.1 <- spGLM(y~X-1, family="poisson", coords=S,
               starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
               tuning=list("beta"=0.1, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
               priors=list("beta.Flat", "phi.Unif"=c(0.03, 10), "sigma.sq.IG"=c(2, 1)),
               amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
               cov.model="exponential", verbose=TRUE, n.report=10)
  t2 <- Sys.time()
  runtime[i, 2] <- t2-t1
  cat("n = ", nseq[i], " time = ", t2-t1, "\n")
}

write.csv(runtime, "spBayes_runtime.txt")