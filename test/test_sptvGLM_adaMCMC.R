rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_sptv_1000.csv")

simdat <- simdat[1:100, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])
time <- as.numeric(simdat$time)

n_postsamp <- 100

mod_out <- sptvGLM_adaMetropGibbs(y = y, X = X, X_tilde = X, 
                                  S = S, time = time, 
                                  family = "poisson", 
                                  N.samp = n_postsamp, 
                                  starting = list(phi_s = c(1, 1), 
                                                  phi_t = c(1, 1),
                                                  beta = c(0, 0)), 
                                  prior = list(phi_s_a = c(0.5, 0.5),
                                               phi_s_b = c(10, 10),
                                               phi_t_a = c(0.5, 0.5),
                                               phi_t_b = c(10, 10),
                                               nu_xi = 1, nu_beta = 3,
                                               nu_z = 3, alpha_epsilon = 0.5),
                                  n.batch = 3, batch.length = 10)


burnin_pc <- 0.5
n_thin <- 5

ids <- 1:n_postsamp
ids <- ids[-(1:(floor(burnin_pc * n_postsamp))+1)]
ids <- ids[c(rep(FALSE, n_thin - 1), TRUE)]

post_beta <- mod_out$beta[, ids]
post_z <- mod_out$z[, ids]
post_phi_s <- mod_out$phi_s[, ids]
post_phi_t <- mod_out$phi_t[, ids]

print(ci_beta(t(post_beta)))

# simdat$postmedian_z <- apply(post_z[1:300, ], 1, median)
# # simdat$postmedian_z2 <- apply(post_z[301:600, ], 1, median)
# leg_title <- TeX('$z(s)$')
# h1 <- 8
# p11 <- pointref_plot(simdat[1:100, ], "z1", legend_title = leg_title, h = h1)
# p12 <- pointref_plot(simdat[100+1:100, ], "z1", legend_title = leg_title, h = h1)
# p13 <- pointref_plot(simdat[200+1:100, ], "z1", legend_title = leg_title, h = h1)
# p21 <- pointref_plot(simdat[1:100, ], "postmedian_z", legend_title = leg_title, h = h1)
# p22 <- pointref_plot(simdat[100+1:100, ], "postmedian_z", legend_title = leg_title, h = h1)
# p23 <- pointref_plot(simdat[200+1:100, ], "postmedian_z", legend_title = leg_title, h = h1)
# gridExtra::grid.arrange(p11, p12, p13, p21, p22, p23, nrow = 2)

post_spParams <- as.data.frame(cbind(t(post_phi_s), t(post_phi_t)))
names(post_spParams) <- c("phi11", "phi12", "phi21", "phi22")
suppressPackageStartupMessages(library(tidyr))
ggplot(gather(post_spParams), aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~key, scales = 'free_x')
