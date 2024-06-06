rm(list = ls())

source("../src/runsrc.R")

# simdat <- read.csv("../data/sim_count1000.csv")
simdat <- read.csv("../data/sim_sptvcount100.3.csv")

# Test on rows 1:100
# simdat <- simdat[1:100, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])
# time <- rep(1, 100)
time <- as.numeric(simdat$time)

n_postsamp <- 500

params <- list(phi = 3.5, phi_t = 0.5, nu_matern = 0.5,
               nu_xi = 1, nu_beta = 2.1, nu_z = 2.1,
               alpha_epsilon = 0.5)

mod_out <- sptGLM_exact(y = y, X = X, S = S, time = time,
                        N.samp = n_postsamp,
                        family = "poisson",
                        mod_params = params)

post_beta <- mod_out$beta
post_z <- mod_out$z

print(ci_beta(post_beta))

simdat$postmedian_z <- apply(post_z[, 1:100], 2, median)
# resids <- tcrossprod(makeG(X), post_z)
# simdat$postmedian_z <- apply(resids, 1, median)

# leg_title <- TeX('$z(s)$')
# p1 <- pointref_plot(simdat, "z1", legend_title = leg_title)
# p6 <- pointref_plot(simdat, "postmedian_z", legend_title = leg_title)
# gridExtra::grid.arrange(p1, p6, ncol = 2)
