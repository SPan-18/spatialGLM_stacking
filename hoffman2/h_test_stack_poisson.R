rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_count1000.csv")

# Test on rows 1:100
simdat <- simdat[1:1000, ]
y <- as.numeric(simdat$y)
X <- as.matrix(simdat[, grep("x", names(simdat))])
S <- as.matrix(simdat[, c("s1", "s2")])

n_postsamp <- 1000

mod_list <- create_model_list(G_decay = c(2.5, 4, 10), 
                              G_smoothness = c(0.5, 1, 1.5),
                              G_epsilon = c(0.25, 0.5),
                              G_nuxi = 0,
                              G_nubeta = 2.1, G_nuz = 2.1)

m_out <- spGLM_stack(y = y, X = X, S = S, N.samp = n_postsamp,
                     family = "poisson",
                     spCov = "matern",
                     solver = "MOSEK",
                     mc.cores = 6,
                     mod_params_list = mod_list)

# elpd_mat <- do.call(cbind, lapply(m_out$models, function(x) x$elpd))
# w_hat <- Mosek_stacking_weights(elpd_mat)

postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
post_z <- postrun_samps$z
post_beta <- postrun_samps$beta
post_xi <- postrun_samps$xi

write.table(post_beta,
            file = "post_stacking/beta.txt",
            col.names = FALSE, row.names = FALSE)
write.table(post_z,
            file = "post_stacking/z.txt",
            col.names = FALSE, row.names = FALSE)
write.table(post_xi,
            file = "post_stacking/xi.txt",
            col.names = FALSE, row.names = FALSE)
