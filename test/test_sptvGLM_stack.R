rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_sptv_1000.csv")

# Test on rows 1:100
simdat <- simdat[1:100, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])
time <- as.numeric(simdat$time)

n_postsamp <- 500

mod_list <- list(G_phi_s = c(2, 3, 5),
                 G_phi_t = c(0.5, 1),
                 G_epsilon = c(0.5),
                 G_nu_xi = c(1),
                 G_nu_beta = c(3),
                 G_nu_z = c(3))

mod_list <- create_candidate_models(mod_list)

set.seed(1729)
m_out <- sptvGLM_stack(y = y, X = X, X_tilde = X, S = S, time = time,
                       N.samp = n_postsamp, MC.samp = 200,
                       family = "poisson",
                       mc.cores = NULL, solver = "MOSEK",
                       mod_params_list = mod_list)

# ci_beta(t(m_out[[1]]$beta))
postrun_samps <- postrunsampler_sptv(m_out, N.samp = n_postsamp)
post_z <- postrun_samps$z
post_xi <- postrun_samps$xi
post_beta <- postrun_samps$beta

print(ci_beta(t(post_beta)))

# simdat$postmedian_z <- apply(post_z[1:100, ], 1, median)
# leg_title <- TeX('$z(s)$')
# h1 <- 8
# p1 <- pointref_plot(simdat, "z1", legend_title = leg_title, h = h1)
# p2 <- pointref_plot(simdat, "postmedian_z", legend_title = leg_title, h = h1)
# gridExtra::grid.arrange(p1, p2, ncol = 2)
