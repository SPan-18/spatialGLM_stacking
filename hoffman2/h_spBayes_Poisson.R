rm(list = ls())

library(spBayes)

# library(ggplot2)
# library(MBA)
# library(latex2exp)
# source("../src/pointrefplot.R")

nseq <- c(1000)
n_run <- length(nseq)

for(i in 1:n_run){
  dat <- read.csv("../data/sim_count1000.csv")
  simdat <- dat[1:nseq[i], ]
  y <- as.numeric(simdat$y)
  X <- as.matrix(simdat[, grep("x", names(simdat))])
  X <- as.matrix(X[, -1])
  S <- as.matrix(simdat[, c("s1", "s2")])
  rm(dat)
  
  ##Collect samples
  beta.starting <- coefficients(glm(y~X, family="poisson"))
  beta.tuning <- t(chol(vcov(glm(y~X, family="poisson"))))
  n.batch <- 500
  batch.length <- 50
  n.samples <- n.batch*batch.length
  
  m.1 <- spGLM(y~X, family="poisson", coords=S,
               starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
               tuning=list("beta"=c(0.1, 0.1), "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
               priors=list("beta.Flat", "phi.Unif"=c(0.03, 10), "sigma.sq.IG"=c(2, 1)),
               amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
               cov.model="exponential", verbose=TRUE, n.report=10)
  
  burn.in <- 0.9*n.samples
  sub.samps <- burn.in:n.samples
  
  w.hat <- m.1$p.w.samples[,sub.samps]
  beta.hat <- m.1$p.beta.theta.samples[sub.samps, 1:2]
  
  write.table(beta.hat, 
              file = paste("post_spBayes/beta_hat_", nseq[i], ".txt", sep = ""), 
              col.names = FALSE, row.names = FALSE)
  write.table(w.hat, 
              file = paste("post_spBayes/z_hat_", nseq[i], ".txt", sep = ""), 
              col.names = FALSE, row.names = FALSE)
}

# w_hat <- read.table("posterior_spatial_Nov28.csv", header = FALSE)
# simdat$postmedian_z <- apply(w_hat, 1, median)
# 
# leg_title <- TeX('$z(s)$')
# p <- pointref_plot(simdat, "postmedian_z", legend_title = leg_title)
# 
# ggsave("postmedian_z_pois.pdf", plot = p, width=4, height=4, units='in')
