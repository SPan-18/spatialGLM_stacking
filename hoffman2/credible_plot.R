# rm(list = ls())

beta_MCMC1 <- read.table("post_MCMC_Mar23/beta.txt")
beta_MCMC1 <- as.data.frame(t(as.matrix(beta_MCMC1)))

beta_stack1 <- read.table("post_stacking/beta.txt")

set.seed(1729)
ids <- sample(1:nrow(beta_MCMC1), 1000)
beta_MCMC1 <- beta_MCMC1[ids, ]
names(beta_MCMC1) = c("beta0", "beta1")

post_beta100 <- beta_MCMC1
post_beta100$Method <- "MCMC"

beta_stack1 <- t(beta_stack1)

rownames(beta_stack1) <- NULL
beta_stack1 <- as.data.frame(beta_stack1)
names(beta_stack1) <- c("beta0", "beta1")
beta_stack1$Method <- "Stacking"

post_beta100 <- rbind.data.frame(post_beta100, beta_stack1, 
                                 stringsAsFactors = TRUE)
row.names(post_beta100) <- NULL
post_beta100$Method <- factor(post_beta100$Method)

# library(reshape2)
library(ggplot2)

p1 <- ggplot(data = post_beta100, aes(x = beta0)) +
  geom_density(aes(fill = Method, color = Method), alpha = 0.3) +
  xlim(0, 10) +
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
        plot.title = element_text(hjust = 0.5, size = 10),
        aspect.ratio = 1)

p4 <- ggplot(MCMCvsStack, aes(x = beta1_MCMC_q, y = beta1_stack_q)) +
  geom_abline(slope = 1, intercept = 0, lty = 3, col = "red") +
  geom_point(col = "grey20", size = 2, alpha = 0.3) +
  labs(title = latex2exp::TeX('$\\beta_1$'),
       x ="MCMC", y = "Stacking") +
  xlim(-0.75, -0.25) +
  ylim(-0.75, -0.25) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        aspect.ratio = 1)

# gridExtra::grid.arrange(p3, p4, ncol = 2)
p34 <- gridExtra::grid.arrange(p3, p4, nrow = 1)

# post_phi <- read.table('post_MCMC_Mar22/phi.txt')
# post_nu <- read.table('post_MCMC_Mar22/nu.txt')
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


# library(ggpubr)
# ggarrange(p34, p_rt, nrow = 1,
#           labels = list("(a)", "(b)"), font.label = list(face = "plain"),
#           common.legend = FALSE, widths = c(0.6, 0.4))

# system2(command = "pdfcrop",
#         args    = c("../plots/credible_runtime.pdf",
#                     "../plots/credible_runtime.pdf")
# )