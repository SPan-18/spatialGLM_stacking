rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_dyncount100.3.csv")

# Test on rows 1:100
# simdat <- simdat[1:1000, ]
y <- as.numeric(simdat$y)
Xt <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])
time <- as.numeric(simdat$time)

n_postsamp <- 500

mod_list <- create_model_list(G_decay = c(3, 4), 
                              G_smoothness = c(0.5, 1),
                              G_timedecay = c(0.5),
                              G_rho = c(0.9),
                              G_epsilon = c(0.5),
                              G_nuxi = 1,
                              G_nubeta = 2.1, G_nuz = 2.1)
set.seed(1729)
m_out <- spDynGLM_stack(y = y, Xt = Xt, S = S, time = time,
                        N.samp = n_postsamp, MC.samp = 200,
                        family = "poisson",
                        spCov = "matern",
                        beta_prior = "nonstationary",
                        omega = 0.95,
                        mc.cores = 4,
                        solver = "MOSEK",
                        mod_params_list = mod_list)

print(ci_beta(t(m_out[[1]][[1]]$beta)))

# postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
# post_z <- postrun_samps$z
# post_xi <- postrun_samps$xi
# post_beta <- postrun_samps$beta

# y.hat <- apply(exp(X %*% post_beta + post_z + post_xi), 2, function(x){rpois(dim(X)[1], x)})
# y.hat.mu <- apply(y.hat, 1, median)

# print(ci_beta(t(post_beta)))
# print(ci_beta(t(m_out$models[[2]]$beta)))

### Compare posterior of z with true z

# simdat$postmedian_z <- apply(post_z, 1, median)
# 
# leg_title <- TeX('$z(s)$')
# title <- TeX('$n = 100$')
# p1 <- pointref_plot(simdat[1:100, ], "z", legend_title = leg_title)
# p2 <- pointref_plot(simdat[1:100,], "postmedian_z", legend_title = leg_title)
# gridExtra::grid.arrange(p1, p2, ncol = 2)
