rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_binom1000.csv")

# Test on rows 1:100
simdat <- simdat[1:100, ]
y <- as.numeric(simdat$y1)
y_binom <- as.numeric(simdat$y2)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])

n_postsamp <- 500
mod_list <- create_model_list(G_decay = c(3, 4), 
                              G_smoothness = 0.5,
                              G_epsilon = 0.5,
                              G_nuxi = 0,
                              G_nubeta = 2.1, G_nuz = 2.1)

m_out <- spGLMM_stack(y = y, X = X, S = S, N.samp = n_postsamp,
                      family = "binomial",
                      n_binom = y_binom,
                      spCov = "matern",
                      mod_params_list = mod_list)

postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
post_z <- postrun_samps$z
post_beta <- postrun_samps$beta

print(ci_beta(t(post_beta)))

simdat$postmedian_z <- apply(post_z, 1, median)
leg_title <- TeX('$z(s)$')
p1 <- pointref_plot(simdat, "z", legend_title = leg_title)
p2 <- pointref_plot(simdat, "postmedian_z", legend_title = leg_title)

ggsave("true_z_binom.pdf", plot = p1, width = 4, height = 4, units='in')
ggsave("postmedian_z_binom.pdf", plot = p2, width = 4, height = 4, units='in')
