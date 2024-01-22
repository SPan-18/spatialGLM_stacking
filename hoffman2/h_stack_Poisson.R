rm(list = ls())

source("../src/runsrc.R")

# n_run <- 14
# n_run <- 2
# nrep <- 5
nseq <- c(100)
# nseq <- c(100, 100)

n_postsamp <- 1000

mod_list <- create_model_list(G_decay = c(3, 5, 10), 
                              G_smoothness = c(0.5, 1, 1.5),
                              G_epsilon = c(0.5, 0.75),
                              G_nuxi = 1,
                              G_nubeta = 2.1, G_nuz = 2.1)

for(i in 1:length(nseq)){
  dat <- read.csv("../data/sim_count1000.csv")
  cat("Running n =", nseq[i], "...\n")
  simdat <- dat[1:nseq[i], ]
  y <- as.numeric(simdat$y)
  X <- as.matrix(simdat[, grep("x", names(simdat))])
  S <- as.matrix(simdat[, c("s1", "s2")])
  rm(dat)
  
  m_out <- spGLM_stack(y = y, X = X, S = S, 
                       N.samp = n_postsamp, MC.samp = 200,
                       family = "poisson",
                       spCov = "matern",
                       mc.cores = 6,
                       solver = "ECOS",
                       mod_params_list = mod_list)
  
  postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
  post_z <- postrun_samps$z
  post_beta <- postrun_samps$beta
  
  write.table(post_beta, 
              file = paste("post_stacking/beta_", nseq[i], ".txt", sep = ""), 
              col.names = FALSE, row.names = FALSE)
  write.table(post_z, 
              file = paste("post_stacking/z_", nseq[i], ".txt", sep = ""), 
              col.names = FALSE, row.names = FALSE)
  
}