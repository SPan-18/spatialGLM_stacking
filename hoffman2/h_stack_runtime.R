rm(list = ls())

source("../src/runsrc.R")

dat <- read.csv("../data/sim_count5000.csv")

n_run <- 14
# n_run <- 2
nrep <- 5
runtime <- array(dim = c(n_run, (nrep+1)))
nseq <- c(c(1:10)*100, 2000, 3000, 4000, 5000)
# nseq <- c(100, 100)
runtime[, 1] <- nseq
colnames(runtime) <- c("n", paste("time", 1:nrep, sep = ""))

n_postsamp <- 500

mod_list <- create_model_list(G_decay = c(3, 4, 10), 
                              G_smoothness = c(0.5, 1, 1.5),
                              G_epsilon = c(0.25, 0.5),
                              G_nuxi = 0,
                              G_nubeta = 2.1, G_nuz = 2.1)

for(i in 1:n_run){
  cat("Running n =", nseq[i], "...\n")
  for(j in 1:nrep){
    simdat <- dat[1:nseq[i], ]
    y <- as.numeric(simdat$y)
    X <- as.matrix(simdat[, grep("x", names(simdat))])
    S <- as.matrix(simdat[, c("s1", "s2")])
    
    t1 <- Sys.time()
    m_out <- try(spGLM_stack(y = y, X = X, S = S, N.samp = n_postsamp,
                         family = "poisson",
                         spCov = "matern",
                         solver = "MOSEK",
                         mc.cores = 6,
                         mod_params_list = mod_list,
                         verbose = FALSE))
    t2 <- Sys.time()
    rt <- difftime(t2, t1, units = "secs")
    runtime[i, (j+1)] <- rt
    write.csv(runtime, "stack_runtime.csv", row.names = FALSE)
    cat("Rep", j, ": Took", rt, units(rt), "\n")
  }
}
