rm(list = ls())

source("../src/h_runsrc.R")

simdat <- read.csv("../data/h_sim_count1000.csv")

# Test on rows 1:100
simdat <- simdat[1:100, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])

n_postsamp <- 500

mod_list <- create_model_list(G_decay = c(3, 4, 10), 
                              G_smoothness = c(0.5, 1, 1.5),
                              G_epsilon = c(0.25, 0.5),
                              G_nuxi = 0,
                              G_nubeta = 2.1, G_nuz = 2.1)

m_out <- spGLM_stack(y = y, X = X, S = S, N.samp = n_postsamp,
                      family = "poisson",
                      spCov = "matern",
                      mod_params_list = mod_list)

postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
post_z <- postrun_samps$z
post_beta <- postrun_samps$beta

# y.hat <- apply(exp(X %*% post_beta + post_z), 2, function(x){rpois(dim(X)[1], x)})
# y.hat.mu <- apply(y.hat, 1, median)

print(ci_beta(t(post_beta)))
# print(ci_beta(t(m_out$models[[2]]$beta)))

### Compare posterior of z with true z

# simdat$postmean_z <- apply(post_z, 1, mean)
# simdat$postsd_z <- apply(post_z, 1, sd)
# simdat$postmedian_z <- apply(post_z, 1, median)
# simdat$yhat <- y.hat.mu
# simdat$postcred1_z <- apply(post_z, 1, function(x) quantile(x, 0.025))
# simdat$postcred2_z <- apply(post_z, 1, function(x) quantile(x, 0.975))
# 
# leg_title <- TeX('$z(s)$')
# p1 <- pointref_plot(simdat, "z", legend_title = leg_title)
# p2 <- pointref_plot(simdat, "postmean_z", legend_title = leg_title)
# p3 <- pointref_plot(simdat, "postsd_z", legend_title = leg_title)
# p4 <- pointref_plot(simdat, "postcred1_z", legend_title = leg_title)
# p5 <- pointref_plot(simdat, "postcred2_z", legend_title = leg_title)
# p6 <- pointref_plot(simdat, "postmedian_z", legend_title = leg_title)
# gridExtra::grid.arrange(p1, p4, p5, ncol = 3)

# pointref_plot(simdat, "y", legend_title = "y")
# pointref_plot(simdat, "yhat", legend_title = "y_hat")

# ggsave("true_z_pois.pdf", plot = p1, width=4, height=4, units='in')
# ggsave("postmedian_z_pois.pdf", plot = p6, width=4, height=4, units='in')

### Histogram of beta
# post_beta <- data.frame(intercept = postrun_samps$beta[1, ],
#                         beta = postrun_samps$beta[2, ])
