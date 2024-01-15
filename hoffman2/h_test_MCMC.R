rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_count1000.csv")

# Test on rows 1:100
simdat <- simdat[1:100, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])

n_postsamp <- 1000000

mod_out <- spGCM_MCMC(y = y, X = X, S = S, 
                      N.samp = n_postsamp,
                      family = "poisson", 
                      spCov = "matern",
                      starting = list(phi = 3, nu = 1,
                                      beta = c(0, 0)),
                      prior = list(phi = c(0.5, 10),
                                   nu = c(0.1, 2),
                                   nu_xi = 1, 
                                   nu_beta = 2.1,
                                   nu_z = 2.1,
                                   alpha_epsilon = 0.5),
                      tuning = c(0.05, 0.05))
                      # n.batch = 10,
                      # batch.length = 10,
                      # n.report = 1)

burnin <- 0.75 * n_postsamp

write.table(mod_out$beta[, -(1:burnin)],
            file = "post_MCMC/beta_100.txt",
            col.names = FALSE, row.names = FALSE)
write.table(mod_out$z[, -(1:burnin)],
            file = "post_MCMC/z_100.txt",
            col.names = FALSE, row.names = FALSE)
write.table(mod_out$xi[, -(1:burnin)],
            file = "post_MCMC/xi_100.txt",
            col.names = FALSE, row.names = FALSE)
write.table(mod_out$phi[-(1:burnin)],
            file = "post_MCMC/phi_100.txt",
            col.names = FALSE, row.names = FALSE)
write.table(mod_out$nu[-(1:burnin)],
            file = "post_MCMC/nu_100.txt",
            col.names = FALSE, row.names = FALSE)


# hist(mod_out$beta[1,])
# hist(mod_out$beta[2,])
# 
# ids <- 1:n_postsamp
# ids <- ids[-(1:floor(0.1 * n_postsamp))]
# ids <- ids[c(rep(FALSE, 9), TRUE)]
# plot(density(mod_out$beta[1, ids]), xlim = c(2, 8))

# plot(density(mod_out$phi))
# plot(density(mod_out$nu))
# 
# plot(mod_out$phi, type = "l")
# plot(mod_out$nu, type = "l")
# plot(mod_out$beta[1,-(1:burnin)], type = "l")
# plot(mod_out$beta[2,], type = "l")

# post_z <- mod_out$z + mod_out$xi
# simdat$postmedian_z <- apply(post_z, 1, median)
# leg_title <- TeX('$z(s)$')
# p1 <- pointref_plot(simdat, "z", legend_title = leg_title)
# p6 <- pointref_plot(simdat, "postmedian_z", legend_title = leg_title)
# gridExtra::grid.arrange(p1, p6, ncol = 2)

