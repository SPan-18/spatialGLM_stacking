rm(list = ls())

source("../src/h_runsrc.R")

simdat <- read.csv("../data/h_sim_count1000.csv")

# Test on rows 1:100
simdat <- simdat[1:200, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])

n_postsamp <- 500
MC_samp <- 200

mod_list <- create_model_list(G_decay = c(3, 4, 10), 
                              G_smoothness = c(0.5, 1, 1.5),
                              G_epsilon = c(0.25, 0.5),
                              G_nuxi = 0,
                              G_nubeta = 2.1, G_nuz = 2.1)

m_out <- spGLM_onlystack(y = y, X = X, S = S, MC.samp = MC_samp,
                         family = "poisson",
                         spCov = "matern",
                         mod_params_list = mod_list)


