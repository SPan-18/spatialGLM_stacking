rm(list = ls())
source("../src/runsrc.R")

# set.seed(1729)
# simdat <- sim_count(n = 100, beta = c(3, 0.5), phi = 3.2)
# write.csv(simdat, "../data/sim_count1000.csv", row.names = FALSE)
simdat <- read.csv("../data/sim_count1000.csv")
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("easting", "northing")])
distmat <- as.matrix(dist(S))

n_postsamp <- 500
mod_list <- create_model_list(G_phi = c(3,4), G_epsilon = 0.75)

t1 <- Sys.time()
m_out <- cm_stacking(y = y, X = X, S = S, N.samp = n_postsamp,
                     mod_params_list = mod_list)
t2 <- Sys.time()
print(t2 - t1)
print(ci_beta(t(m_out[[1]]$post_samples$beta)))
print(ci_beta(t(m_out[[2]]$post_samples$beta)))


### Compare plot of z!!!

