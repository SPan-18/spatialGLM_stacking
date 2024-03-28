rm(list = ls())

source("../src/runsrc.R")

simdat <- read.csv("../data/sim_count1000.csv")

post_z_df <- read.table("post_MCMC_Mar23/z.txt")
post_xi_df <- read.table("post_MCMC_Mar23/xi.txt")
post_z <- as.matrix(post_z_df)
post_xi <- as.matrix(post_xi_df)
post_sp <- post_z + post_xi
simdat$postmedian_z_mcmc <- apply(post_sp, 1, median)

post_z_df <- read.table("post_stacking/z.txt")
post_xi_df <- read.table("post_stacking/xi.txt")
post_z <- as.matrix(post_z_df)
post_xi <- as.matrix(post_xi_df)
post_sp <- post_z + post_xi
simdat$postmedian_z_stack <- apply(post_sp, 1, median)

p1 <- pointref_plot(simdat, "z", legend_title = NULL)
p2 <- pointref_plot(simdat, "postmedian_z_mcmc", legend_title = NULL)
p3 <- pointref_plot(simdat, "postmedian_z_stack", legend_title = NULL)

library(ggpubr)
ggarrange(p1, p2, p3, nrow = 1, 
          labels = list("(a)", "(b)", "(c)"), font.label = list(face = "plain"),
          common.legend = TRUE, legend="right")

# # after saving image, run pdfcrop
# system2(command = "pdfcrop",
#         args    = c("../plots/surface_plots.pdf",
#                     "../plots/surface_plots.pdf")
# )
