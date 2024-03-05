rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_dyncount100.10.csv")

# Test on rows 1:100
# simdat <- simdat[1:1000, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])
time <- as.numeric(simdat$time)

n_postsamp <- 500

params <- list(phi = 3.5, phi_t = 0.5, nu_matern = 0.5, rho = 0.5, 
               nu_xi = 1, nu_beta = 2.1, nu_z = 2.1,
               alpha_epsilon = 0.5)

mod_out <- spDynGLM_exact(y = y, X = X, S = S, time = time,
                          N.samp = n_postsamp,
                          family = "poisson",
                          omega = 0.95,
                          mod_params = params)

post_beta <- mod_out$beta
post_z <- mod_out$z

print(ci_beta(post_beta))


# simdat$postmedian_z <- apply(post_z, 2, median)
# leg_title <- TeX('$z(s)$')
# p1 <- pointref_plot(simdat, "z", legend_title = leg_title)
# p2 <- pointref_plot(simdat, "postmedian_z", legend_title = leg_title)
# gridExtra::grid.arrange(p1, p2, ncol = 2)

summary_beta <- ci_beta(post_beta)
beta_true <- c(0.5, 0.7, 0.9, 1.2, 1.7, 2.2, 3, 4, 5.5, 7, 
               -1.5, 0.2, 1.3, -0.3, 0.01, 0.2, -0.5, 0.12, -0.6, 0.26)

ids <- 1:10
ymin1 <- min(summary_beta[ids, '2.5%'])
ymax1 <- max(summary_beta[ids, '97.5%'])

plot(1:length(ids), summary_beta[ids, '50%'], col = "steelblue4", lwd = 2,
     type = "l", ylim = c(ymin1, ymax1),
     xlab = latex2exp::TeX('$t$'), 
     ylab = latex2exp::TeX('$\\beta_1(t)$'))
polygon(c(rev(1:length(ids)), 1:length(ids)),
        c(rev(summary_beta[ids, "97.5%"]),
          summary_beta[ids, "2.5%"]),
        lty = 2,
        border = NA,
        col = adjustcolor("steelblue1", alpha.f = 0.2))
lines(1:length(ids), summary_beta[ids, "2.5%"], lwd = 1, col = "steelblue1")
lines(1:length(ids), summary_beta[ids, "97.5%"], lwd = 1, col = "steelblue1")
lines(1:length(ids), beta_true[ids], lwd = 2, col = "darkred", lty = 2)
legend(7.5, 2.2, legend = c("Truth", "Posterior"),
       col = c("darkred", "steelblue1"), lty = c(2, 1), lwd = c(2, 2),
       bty = "n", y.intersp = 0.25)

