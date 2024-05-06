rm(list = ls())
library(tidyverse)

source("../src/runsrc.R")

bbs21 <- read.csv("../data/bbs/BBS21.csv")
bbs21 <- bbs21 %>%
  filter(Latitude < 60) %>%
  filter(!(RouteDataID == 6376050)) %>%
  filter(BirdCount < 500)

ids <- 1:nrow(bbs21)
y <- as.numeric(bbs21[ids, "BirdCount"])
X <- as.matrix(bbs21[, c("Longitude", "Latitude", "NCar", "Noise")])
X <- cbind(rep(1, length(y)), X)
S <- as.matrix(bbs21[ids, c("Longitude", "Latitude")])

n_postsamp <- 20000

mod_out <- spGCM_adaMetropGibbs(y = y, X = X, S = S, 
                                family = "poisson", 
                                N.samp = n_postsamp,
                                spCov = "matern",
                                starting = list(phi = 500, nu = 1,
                                                beta = rep(0, 5)),
                                prior = list(phi = c(0.5, 1500),
                                             nu = c(0.1, 2),
                                             nu_xi = 1, 
                                             nu_beta = 2.1,
                                             nu_z = 2.1,
                                             alpha_epsilon = 0.5))

ids <- 1:n_postsamp
# ids <- ids[-(1:(floor(0.1 * n_postsamp))+1)]
# ids <- ids[c(rep(FALSE, 8), TRUE)]
# 
write.table(mod_out$beta[, ids],
            file = "bbs21_MCMC/beta.txt",
            col.names = FALSE, row.names = FALSE)
write.table(mod_out$z[, ids],
            file = "bbs21_MCMC/z.txt",
            col.names = FALSE, row.names = FALSE)
write.table(mod_out$xi[, ids],
            file = "bbs21_MCMC/xi.txt",
            col.names = FALSE, row.names = FALSE)
write.table(mod_out$phi[ids],
            file = "bbs21_MCMC/phi.txt",
            col.names = FALSE, row.names = FALSE)
write.table(mod_out$nu[ids],
            file = "bbs21_MCMC/nu.txt",
            col.names = FALSE, row.names = FALSE)
