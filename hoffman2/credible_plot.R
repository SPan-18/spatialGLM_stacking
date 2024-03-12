rm(list = ls())

# beta_MCMC1 <- read.table("post_spBayes/beta_hat_100.txt")
beta_MCMC1 <- read.table("post_MCMC_jan21/beta_100.txt")
beta_MCMC1 <- as.data.frame(t(as.matrix(beta_MCMC1)))
# beta_MCMC2 <- read.table("post_spBayes/beta_hat_500.txt")

beta_stack1 <- read.table("post_stacking/beta_100.txt")
# beta_stack2 <- read.table("post_stacking/beta_500.txt")
# beta_stack3 <- read.table("post_stacking/beta_1000.txt")

ids <- sample(1:nrow(beta_MCMC1), 1000)
beta_MCMC1 <- beta_MCMC1[ids, ]
# beta_MCMC2 <- beta_MCMC2[ids, ]
names(beta_MCMC1) = c("beta0", "beta1")
# names(beta_MCMC2) = c("beta0", "beta1")

post_beta100 <- beta_MCMC1
# post_beta500 <- beta_MCMC2
# post_beta1000 <- beta_MCMC3
post_beta100$Method <- "MCMC"
# post_beta500$Method <- "MCMC"
# post_beta1000$Method <- "MCMC"

beta_stack1 <- t(beta_stack1)
# beta_stack2 <- t(beta_stack2)
# beta_stack3 <- t(beta_stack3)

rownames(beta_stack1) <- NULL
# rownames(beta_stack2) <- NULL
# rownames(beta_stack3) <- NULL
beta_stack1 <- as.data.frame(beta_stack1)
# beta_stack2 <- as.data.frame(beta_stack2)
# beta_stack3 <- as.data.frame(beta_stack3)
names(beta_stack1) <- c("beta0", "beta1")
# names(beta_stack2) <- c("beta0", "beta1")
# names(beta_stack3) <- c("beta0", "beta1")
beta_stack1$Method <- "Stacking"
# beta_stack2$Method <- "Stacking"
# beta_stack3$Method <- "Stacking"

post_beta100 <- rbind.data.frame(post_beta100, beta_stack1, 
                                 stringsAsFactors = TRUE)
# post_beta500 <- rbind.data.frame(post_beta500, beta_stack2, 
#                                  stringsAsFactors = TRUE)
# post_beta1000 <- rbind.data.frame(post_beta1000, beta_stack3, 
#                                  stringsAsFactors = TRUE)
row.names(post_beta100) <- NULL
# row.names(post_beta500) <- NULL
# row.names(post_beta1000) <- NULL
post_beta100$Method <- factor(post_beta100$Method)
# post_beta500$Method <- factor(post_beta500$Method)
# post_beta1000$Method <- factor(post_beta1000$Method)

# library(reshape2)
library(ggplot2)

# post_beta100 <- melt(post_beta100, id.vars = "Method")
# 
# post_beta100 <- data.table::data.table(post_beta100)
# post_beta100[variable == "beta0", x_min := 0]
# post_beta100[variable == "beta0", x_max := 10]

# pal <- wesanderson::wes_palette("Darjeeling1", 2, type = "discrete")
p1 <- ggplot(data = post_beta100, aes(x = beta0)) +
  geom_density(aes(fill = Method, color = Method), alpha = 0.3) +
  xlim(0, 10) +
  # scale_color_manual(values = pal) +
  # scale_fill_manual(values = pal) +
  theme_bw() +
  xlab(latex2exp::TeX('$\\beta_0$')) +
    geom_point(aes(x = 5, y = 0),
               pch = 24, fill = "#B4464B", size = 2,
               inherit.aes = F) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.83, 0.9),
        legend.background = element_blank(),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size = 8, face = "italic"),
        legend.title.align = 0.5,
        legend.text = element_text(size = 9),
        axis.title.y=element_blank(),
        aspect.ratio = 1)

p2 <- ggplot(data = post_beta100, aes(x = beta1)) +
  geom_density(aes(fill = Method, color = Method), alpha = 0.3) +
  # scale_color_manual(values = c("#B4AF46", "#4682B4")) +
  # scale_fill_manual(values = c("#B4AF46", "#4682B4")) +
  xlim(-2.5, 1.5) +
  theme_bw() +
  xlab(latex2exp::TeX('$\\beta_1$')) +
  geom_point(aes(x = -0.5, y = 0),
             pch = 24, fill = "#B4464B", size = 2,
             inherit.aes = F) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.85, 0.9),
        legend.background = element_blank(),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size = 8, face = "italic"),
        legend.title.align = 0.5,
        legend.text = element_text(size = 9),
        axis.title.y=element_blank(),
        aspect.ratio = 1)

# gridExtra::grid.arrange(p1, p2, ncol = 2)


# post_beta100 <- reshape2::melt(post_beta100, id.vars = "Method")
# p2 <- ggplot(data = post_beta100, aes(x = value)) +
#   geom_density(aes(fill = Method), alpha = 0.2) +
#   facet_wrap(~ variable, scales = "free_x")

# qq <- as.data.frame(qqplot(beta_MCMC1[, 1], 
#                            beta_stack1$beta0, 
#                            plot.it = FALSE))
# ggplot(qq) + 
#   geom_point(aes(x = x, y = y))



# phi_MCMC <- read.table("post_MCMC_jan21/phi_100.txt")
# nu_MCMC <- read.table("post_MCMC_jan21/nu_100.txt")
# 
# plot(density(phi_MCMC[, 1]))
# plot(density(nu_MCMC[, 1]))


beta0_MCMC <- beta_MCMC1[, 1]
beta1_MCMC <- beta_MCMC1[, 2]

beta0_stack <- beta_stack1$beta0
beta1_stack <- beta_stack1$beta1

qseq <- seq(0, 1, length.out = 100)
beta0_MCMC_q <- quantile(beta0_MCMC, qseq)
beta0_stack_q <- quantile(beta0_stack, qseq)
beta1_MCMC_q <- quantile(beta1_MCMC, qseq)
beta1_stack_q <- quantile(beta1_stack, qseq)

MCMCvsStack <- data.frame(beta0_MCMC_q = beta0_MCMC_q,
                          beta0_stack_q = beta0_stack_q,
                          beta1_MCMC_q = beta1_MCMC_q,
                          beta1_stack_q = beta1_stack_q)

p3 <- ggplot(MCMCvsStack, aes(x = beta0_MCMC_q, y = beta0_stack_q)) +
  geom_abline(slope = 1, intercept = 0, lty = 3, col = "red") +
  geom_point(col = "grey20", size = 2, alpha = 0.3) +
  labs(title="(Intercept)",
       x ="MCMC", y = "Stacking") +
  xlim(0, 10) +
  ylim(0, 10) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1)

p4 <- ggplot(MCMCvsStack, aes(x = beta1_MCMC_q, y = beta1_stack_q)) +
  geom_abline(slope = 1, intercept = 0, lty = 3, col = "red") +
  geom_point(col = "grey20", size = 2, alpha = 0.3) +
  labs(title = latex2exp::TeX('$\\beta_1$'),
       x ="MCMC", y = "Stacking") +
  xlim(-1.5, 0.5) +
  ylim(-1.5, 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1)

# gridExtra::grid.arrange(p3, p4, ncol = 2)
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 4)

# plot(beta0_MCMC_q, beta0_stack_q, xlim = c(0, 10), ylim = c(0, 10))
# abline(0, 1)
# 
# plot(beta1_MCMC_q, beta1_stack_q)
# abline(0, 1)

# post_phi <- read.table('post_MCMC_jan21/phi_100.txt')
# post_nu <- read.table('post_MCMC_jan21/nu_100.txt')
# 
# LaplacesDemon::MCSE(beta1_MCMC)
# LaplacesDemon::MCSE(beta0_MCMC)
# LaplacesDemon::MCSE(post_phi[, 1])
# LaplacesDemon::MCSE(post_nu[, 1])
# 
# LaplacesDemon::ESS(beta1_MCMC)
# LaplacesDemon::ESS(beta0_MCMC)
# LaplacesDemon::ESS(post_phi[, 1])
# LaplacesDemon::ESS(post_nu[, 1])



