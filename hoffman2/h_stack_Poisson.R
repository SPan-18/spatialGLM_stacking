rm(list = ls())

source("../src/runsrc.R")

n_run <- 14
# n_run <- 2
nrep <- 5
nseq <- c(100, 200, 500, 1000)
# nseq <- c(100, 100)

n_postsamp <- 1000

mod_list <- create_model_list(G_decay = c(3, 4, 10), 
                              G_smoothness = c(0.5, 1, 1.5),
                              G_epsilon = c(0.25, 0.5),
                              G_nuxi = 0,
                              G_nubeta = 2.1, G_nuz = 2.1)

for(i in 1:n_run){
  dat <- read.csv("../data/sim_count1000.csv")
  cat("Running n =", nseq[i], "...\n")
  simdat <- dat[1:nseq[i], ]
  y <- as.numeric(simdat$y)
  X <- as.matrix(simdat[, grep("x", names(simdat))])
  S <- as.matrix(simdat[, c("s1", "s2")])
  rm(dat)
  
  m_out <- try(spGLM_stack(y = y, X = X, S = S, N.samp = n_postsamp,
                           family = "poisson",
                           spCov = "matern",
                           solver = "MOSEK",
                           mc.cores = 6,
                           mod_params_list = mod_list,
                           verbose = FALSE))
  
  #### write code for finding MLPD for hold-out
  
}