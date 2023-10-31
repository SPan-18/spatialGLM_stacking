rm(list = ls())
source("../src/runsrc.R")

simdat <- sim_count(n = 100, beta = c(3, 0.5), phi = 3.2)
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("easting", "northing")])
distmat <- as.matrix(dist(S))

n_postsamp <- 100
mod_list <- create_model_list(G_phi = c(3, 4), G_epsilon = 0.75)

# m_out <- CV_posterior_sampler(y = y, X = X, distmat = distmat, N.samp = n_postsamp,
#                               family = "poisson",
#                               beta_prior = "gaussian",
#                               spatial_prior = "gaussian",
#                               mod_params = mod_list[[1]],
#                               CV_K = 10,
#                               Rfast_parallel = FALSE)
# ci_beta(m_out[[2]]$beta)

m_out <- cm_stacking(y = y, X = X, S = S, N.samp = n_postsamp,
                     mod_params_list = mod_list)
