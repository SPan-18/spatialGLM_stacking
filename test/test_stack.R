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
mod_list <- create_model_list(G_phi = c(2, 4.5), G_epsilon = c(0.5, 0.75))

m_out <- cm_stacking(y = y, X = X, S = S, N.samp = n_postsamp,
                     mod_params_list = mod_list)

post_z <- postrunsampler_z(m_out, 1000)

# print(ci_beta(t(m_out$models[[1]]$beta)))
# print(ci_beta(t(m_out[[2]]$post_samples$beta)))

### Compare posterior of z with true z

# simdat$postmean_z <- apply(post_z, 1, mean)
# simdat$postsd_z <- apply(post_z, 1, sd)
# 
# leg_title <- TeX('$z(s)$')
# p1 <- pointref_plot(simdat, "z", legend_title = leg_title)
# p2 <- pointref_plot(simdat, "postmean_z", legend_title = leg_title)
# p3 <- pointref_plot(simdat, "postsd_z", legend_title = leg_title)
# gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
