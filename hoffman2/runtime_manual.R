library(latex2exp)
library(ggplot2)

run_time_stack <- read.csv("stack_runtime_Nov16.csv")
run_time_MCMC <- read.csv("spBayes_runtime_Nov13.csv")
nseq <- run_time$n

stack_fit <- c(15.1, 25.2, 38.3, 57.2, 75.6) # 99, 115.8, 138.0, 156.0, 198
# stack_nofit <- c(10.3, 17.8, 25.6, 40.8, 45.1) # 55.7, 73.8, 90.6, 105.6, 126.6
spBayes_time <- c(55.4, 381.6, 981.6, 2469, 4968)

run_time <- data.frame(n = nseq, stack_nofit = stack_nofit,
                       stack_fit = stack_fit, spBayes_time = spBayes_time)

library("tidyverse")
rt <- run_time %>%
  select(n, stack_nofit, stack_fit, spBayes_time) %>%
  gather(key = "variable", value = "value", -n)

ggplot(rt, aes(x = n, y = value)) +
  geom_line(aes(color = variable)) +
  geom_point(aes(color = variable, alpha = 0.5)) +
  scale_y_continuous(trans='log10') +
  xlab(TeX('$n$')) +
  ylab(TeX('$\\log_{10}$time')) +
  labs(color = "Method") +
  scale_color_discrete(labels = c("spBayes", "stack.fit", "stack")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.12, 0.99),
        legend.text = element_text(size = 10, family = "mono"),
        legend.key.size = unit(0.5, 'cm'),
        axis.title.x = element_text(size = 12, family = "mono"),
        axis.title.y = element_text(size = 12, family = "mono"),
        # axis.line = element_line(color='black'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        aspect.ratio = 0.67)
