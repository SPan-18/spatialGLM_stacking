rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_sptmpcount100.3.csv")

# Test on rows 1:100
# simdat <- simdat[1:1000, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])
time <- as.numeric(simdat$time)

n_postsamp <- 500

mod_list <- create_model_list(G_decay = c(3, 4), 
                              G_smoothness = c(0.5, 1),
                              G_timedecay = c(0.5),
                              G_epsilon = c(0.5),
                              G_nuxi = 1,
                              G_nubeta = 2.1, G_nuz = 2.1)
set.seed(1729)
m_out <- sptGLM_stack(y = y, X = X, S = S, time = time,
                        N.samp = n_postsamp, MC.samp = 200,
                        family = "poisson",
                        spCov = "matern",
                        mc.cores = 4,
                        solver = "MOSEK",
                        mod_params_list = mod_list)

# print(ci_beta(t(m_out[[1]][[1]]$beta)))
postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
post_z <- postrun_samps$z
post_xi <- postrun_samps$xi
post_beta <- postrun_samps$beta

print(ci_beta(t(post_beta)))

# simdat$postmedian_z <- apply(post_z + post_xi, 1, median)
# leg_title <- TeX('$z(s)$')
# p1 <- pointref_plot(simdat, "z", legend_title = leg_title)
# p2 <- pointref_plot(simdat, "postmedian_z", legend_title = leg_title)
# gridExtra::grid.arrange(p1, p2, ncol = 2)
