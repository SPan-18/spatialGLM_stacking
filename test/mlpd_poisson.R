rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_count1000.csv")

# sample size
N <- 100
# Holdout rows 1:n_h
n_h <- 100
holdout <- simdat[1:n_h, ]
simdat <- simdat[n_h + 1:N, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])

n_postsamp <- 500

mod_list <- create_model_list(G_decay = c(3, 4, 10), 
                              G_smoothness = c(0.5, 1, 1.5),
                              G_epsilon = c(0.25, 0.5),
                              G_nuxi = 0,
                              G_nubeta = 2.1, G_nuz = 2.1)

m_out <- spGLM_stack(y = y, X = X, S = S, N.samp = n_postsamp,
                     family = "poisson",
                     spCov = "matern",
                     mc.cores = 6,
                     solver = "MOSEK",
                     mod_params_list = mod_list)

postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
post_z <- postrun_samps$z
post_beta <- postrun_samps$beta

