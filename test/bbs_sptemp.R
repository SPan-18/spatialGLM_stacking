rm(list = ls())

source("../src/runsrc.R")

bbs <- read.csv("../data/BBS10-19.csv")
bbs <- bbs %>%
  filter(Latitude < 55) %>%
  filter(BirdCount < 500)

# Test on rows 1:100
# simdat <- simdat[1:1000, ]
y <- as.numeric(bbs[, "BirdCount"])
Xt <- as.matrix(bbs[, c("NCar", "Noise")])
Xt <- cbind(rep(1, length(y)), Xt)
S <- as.matrix(bbs[, c("Longitude", "Latitude")])
time <- as.numeric(bbs$Time)

n_postsamp <- 100

mod_list <- create_model_list(G_decay = c(400, 800), 
                              G_smoothness = c(0.5, 1),
                              G_timedecay = c(0.9),
                              G_rho = c(0.9),
                              G_epsilon = c(0.5),
                              G_nuxi = 1,
                              G_nubeta = 2.1, G_nuz = 2.1)
set.seed(1729)
m_out <- spDynGLM_stack(y = y, Xt = Xt, S = S, time = time,
                        N.samp = n_postsamp, MC.samp = 200,
                        family = "poisson",
                        spCov = "matern",
                        beta_prior = "stationary",
                        omega = NULL,
                        mc.cores = 4,
                        solver = "MOSEK",
                        mod_params_list = mod_list)

postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
post_z <- postrun_samps$z
post_xi <- postrun_samps$xi
post_beta <- postrun_samps$beta

summary_beta <- ci_beta(t(post_beta))
# print(summary_beta)

ids <- 1:10 + 10 * 2
ymin1 <- min(summary_beta[ids, '2.5%'])
ymax1 <- max(summary_beta[ids, '97.5%'])

xtime <- 2010:2019
plot(xtime, summary_beta[ids, '50%'], col = "steelblue4", lwd = 2,
     type = "l", ylim = c(ymin1, ymax1),
     xlab = latex2exp::TeX('Year'), 
     ylab = "Noise")
polygon(c(rev(xtime), xtime),
        c(rev(summary_beta[ids, "97.5%"]),
          summary_beta[ids, "2.5%"]),
        lty = 2,
        border = NA,
        col = adjustcolor("steelblue1", alpha.f = 0.2))
# lines(1:length(ids), summary_beta[ids, "2.5%"], lwd = 1, col = "steelblue1")
# lines(1:length(ids), summary_beta[ids, "97.5%"], lwd = 1, col = "steelblue1")
