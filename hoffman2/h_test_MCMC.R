rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_count1000.csv")

# Test on rows 1:100
simdat <- simdat[1:100, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])

n_postsamp <- 1000

mod_out <- spGCM_adaMetropGibbs(y = y, X = X, S = S, 
                                family = "poisson", 
                                N.samp = n_postsamp,
                                spCov = "matern",
                                starting = list(phi = 3, nu = 1,
                                                beta = c(0, 0)),
                                prior = list(phi = c(0.5, 10),
                                             nu = c(0.1, 2),
                                             nu_xi = 1, 
                                             nu_beta = 2.1,
                                             nu_z = 2.1,
                                             alpha_epsilon = 0.5))


ids <- 1:n_postsamp
# ids <- ids[-(1:(floor(0.1 * n_postsamp))+1)]
# ids <- ids[c(rep(FALSE, 8), TRUE)]

# 
# write.table(mod_out$beta[, ids],
#             file = "post_MCMC/beta.txt",
#             col.names = FALSE, row.names = FALSE)
# write.table(mod_out$z[, ids],
#             file = "post_MCMC/z.txt",
#             col.names = FALSE, row.names = FALSE)
# write.table(mod_out$xi[, ids],
#             file = "post_MCMC/xi.txt",
#             col.names = FALSE, row.names = FALSE)
# write.table(mod_out$phi[ids],
#             file = "post_MCMC/phi.txt",
#             col.names = FALSE, row.names = FALSE)
# write.table(mod_out$nu[ids],
#             file = "post_MCMC/nu.txt",
#             col.names = FALSE, row.names = FALSE)
