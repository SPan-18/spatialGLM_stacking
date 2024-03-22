rm(list = ls())

library(latex2exp)
library(ggplot2)

run_time_stack <- read.csv("stack_runtime_Nov28.csv")
run_time_MCMC <- read.csv("spBayes_runtime_Nov27.csv")
nseq <- run_time_stack$n

# stack_fit <- c(15.1, 25.2, 38.3, 57.2, 75.6) # 99, 115.8, 138.0, 156.0, 198
# stack_nofit <- c(10.3, 17.8, 25.6, 40.8, 45.1) # 55.7, 73.8, 90.6, 105.6, 126.6
# spBayes_time <- c(55.4, 381.6, 981.6, 2469, 4968)   

stack_fit <- apply(as.matrix(run_time_stack[, -1]), 1, mean)
spBayes_time <- run_time_MCMC$time

# subset.ids <- 1:14
# run_time <- data.frame(n = nseq[subset.ids], # stack_nofit = stack_nofit,
#                        stack_fit = stack_fit[subset.ids], 
#                        spBayes_time = c(spBayes_time, rep(max(spBayes_time), 4))[subset.ids])
# 

nseq_pred <- c(nseq, seq(110, 4990, length.out = 50))
y1 <- stack_fit
y2 <- spBayes_time
x1 <- nseq
x2 <- nseq[1:10]
mod_stack <- lm(y1 ~ poly(x1, 3, raw = TRUE))
stack_rt_pred <- predict(mod_stack, data.frame(x1 = nseq_pred))

mod_spBayes <- lm(y2 ~ poly(x2, 3, raw = TRUE))
spBayes_rt_pred <- predict(mod_spBayes, data.frame(x2 = nseq_pred))

run_time_obs <- data.frame(n = nseq, stack_rt = stack_rt_pred[1:length(nseq)],
                           spBayes_rt = c(spBayes_rt_pred[1:10], rep(NA, 4)))
run_time_pred <- data.frame(n = nseq_pred, stack_rt = stack_rt_pred, 
                            spBayes_rt = spBayes_rt_pred)

library("tidyverse")
rt_obs <- run_time_obs %>%
  select(n, stack_rt, spBayes_rt) %>%
  gather(key = "variable", value = "value", -n) %>%
  na.omit()
rt_df <- run_time_pred %>%
  select(n, stack_rt, spBayes_rt) %>%
  gather(key = "variable", value = "value", -n)

ggplot(rt_df, aes(x = n, linetype = variable)) +
  geom_line(aes(y = value)) +
  # geom_point(data = rt_obs, size = 1, aes(x = n, y = value)) +
  # scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans='log10') +
  xlab('Sample size') +
  ylab(TeX('$\\log_{10}$Time')) +
  # labs(linetype = "Method") +
  # scale_color_discrete(labels = c("MCMC", "Stacking")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme_bw() +
  theme(legend.position="none",
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        # axis.line = element_line(color='black'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        aspect.ratio = 0.66667)

# ggplot(rt_df, aes(x = n, color = variable)) +
#   geom_line(aes(y = value)) +
#   geom_point(data = rt_obs, aes(x = n, y = value)) +
#   # scale_x_continuous(trans = 'log10') +
#   scale_y_continuous(trans='log10') +
#   xlab('Sample size') +
#   ylab(TeX('$\\log_{10}$Time')) +
#   labs(color = "Method") +
#   scale_color_discrete(labels = c("MCMC", "Stacking")) +
#   theme_bw() +
#   theme(legend.title = element_text(size = 9, face = "italic"),
#         legend.title.align=0.5,
#         legend.background = element_blank(),
#         legend.position = c(0.12, 0.9),
#         legend.text = element_text(size = 10),
#         legend.key.size = unit(0.3, 'cm'),
#         axis.title.x = element_text(size = 11),
#         axis.title.y = element_text(size = 11),
#         # axis.line = element_line(color='black'),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         aspect.ratio = 0.66667)
