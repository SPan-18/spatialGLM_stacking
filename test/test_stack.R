rm(list = ls())
source("../src/runsrc.R")

# set.seed(1729)
# simdat <- sim_count(n = 1000, beta = c(3, 0.5), phi = 3.2)
# write.csv(simdat, "../data/sim_count1000.csv", row.names = FALSE)
simdat <- read.csv("../data/sim_count1000.csv")
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])
distmat <- as.matrix(dist(S))

n_postsamp <- 500
mod_list <- create_model_list(G_phi = c(2, 3.5), G_epsilon = c(0.5, 0.75))

m_out <- cm_stacking(y = y, X = X, S = S, N.samp = n_postsamp,
                     mod_params_list = mod_list)

postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
post_z <- postrun_samps$z

# print(ci_beta(t(m_out$models[[1]]$beta)))
# print(ci_beta(t(m_out$models[[2]]$beta)))

### Compare posterior of z with true z

# simdat$postmean_z <- apply(post_z, 1, mean)
# simdat$postsd_z <- apply(post_z, 1, sd)
simdat$postmedian_z <- apply(post_z, 1, median)
# simdat$postcred1_z <- apply(post_z, 1, function(x) quantile(x, 0.025))
# simdat$postcred2_z <- apply(post_z, 1, function(x) quantile(x, 0.975))
# 
leg_title <- TeX('$z(s)$')
# p1 <- pointref_plot(simdat, "z", legend_title = leg_title)
# p2 <- pointref_plot(simdat, "postmean_z", legend_title = leg_title)
# p3 <- pointref_plot(simdat, "postsd_z", legend_title = leg_title)
# p4 <- pointref_plot(simdat, "postcred1_z", legend_title = leg_title)
# p5 <- pointref_plot(simdat, "postcred2_z", legend_title = leg_title)
p6 <- pointref_plot(simdat, "postmedian_z", legend_title = leg_title)
# gridExtra::grid.arrange(p1, p4, p5, ncol = 3)

### Histogram of beta
# post_beta <- data.frame(intercept = postrun_samps$beta[1, ],
#                         beta = postrun_samps$beta[2, ])
