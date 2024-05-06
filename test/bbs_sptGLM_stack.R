rm(list = ls())

source("../src/runsrc.R")

bbs <- read.csv("../data/BBS10-19.csv")
bbs <- bbs %>%
  filter(Latitude < 55) %>%
  filter(BirdCount < 500)

# Test on rows 1:100
# simdat <- simdat[1:1000, ]
y <- as.numeric(bbs[, "BirdCount"])
X <- as.matrix(bbs[, c("NCar", "Noise")])
X <- cbind(rep(1, length(y)), X)
S <- as.matrix(bbs[, c("Longitude", "Latitude")])
time <- as.numeric(bbs$Time)

n_postsamp <- 500

mod_list <- create_model_list(G_decay = c(400, 800),
                              G_smoothness = c(0.5, 1),
                              G_timedecay = c(0.5),
                              G_epsilon = c(0.5),
                              G_nuxi = 1,
                              G_nubeta = 3, G_nuz = 3)

set.seed(1729)
m_out <- sptGLM_stack(y = y, X = X, S = S, time = time,
                      N.samp = n_postsamp, MC.samp = 200,
                      family = "poisson",
                      spCov = "matern",
                      mc.cores = NULL,
                      solver = "MOSEK",
                      mod_params_list = mod_list)

postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
post_z <- postrun_samps$z
post_xi <- postrun_samps$xi
post_beta <- postrun_samps$beta

print(ci_beta(t(post_beta)))

# bbs$postmedian_z <- apply(post_z, 1, median)
# 
# library("sf")
# library("rnaturalearth")
# library("rnaturalearthdata")
# library("sp")
# library("raster")
# library("ggpubr")
# 
# usa <- ne_countries(country = "United States of America",
#                     scale = "medium", returnclass = "sf")
# # usa_map <- map('usa')
# # usa <- st_as_sf(usa_map)
# latmin <- min(bbs$Latitude)
# latmax <- max(bbs$Latitude)
# lonmin <- min(bbs$Longitude)
# lonmax <- max(bbs$Longitude)
# lat_extra <- abs(latmax - latmin) * 0.00
# lon_extra <- abs(lonmax - lonmin) * 0.00
# latlim <- c(24.5977, 49.6)
# lonlim <- c(-125, -66.5)
# 
# blue.red <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# col.br <- colorRampPalette(blue.red)
# brks <- c(-1.0, -0.5, 0, 0.5, 1.0)
# 
# plot_list <- vector(mode = "list", length = 10)
# 
# for(i in 1:10){
#   ids <- which(bbs$Time == i)
#   surf <- mba.surf(bbs[ids , c("Longitude","Latitude", "postmedian_z")],
#                    no.X = 200, no.Y = 200, h = 6, m = 1, n = 1,
#                    extend=FALSE)$xyz.est
#   surf_rast <- raster(surf)
#   us_rast <- raster::extract(surf_rast, usa, cellnumbers = TRUE)
#   full_grid <- expand.grid(surf$x, rev(surf$y))
#   us_df <- data.frame(full_grid[us_rast[[1]][, "cell"], ],
#                       z = us_rast[[1]][, "value"])
#   colnames(us_df) <- c("x", "y", "z")
#   us_df <- na.omit(us_df)
#   
#   plot_list[[i]] <- ggplot(usa) +
#     geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
#     scale_fill_gradientn(colours = col.br(20), breaks = brks, 
#                          label = function(x) sprintf("%.1f", x)) +
#     geom_sf(data = usa, color = "black", fill = NA, linewidth = 0.05) +
#     coord_sf(xlim = lonlim, ylim = latlim, expand = FALSE) +
#     labs(fill = TeX('$z(s, t)$'), 
#          title = as.character(2009+i)) + 
#     theme_bw() + 
#     theme(legend.position = "none",
#           panel.grid = element_blank(), 
#           panel.border = element_blank(),
#           panel.background = element_blank(), 
#           axis.title=element_blank(),
#           axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           plot.title = element_text(color = "gray40", size = 10,
#                                     hjust = 0.5, vjust = -1.5),
#           plot.margin = unit(c(0,0,0.25,0), "cm"),
#           aspect.ratio = 0.75)
# }
# 
# p1 <- ggarrange(plot_list[[1]], plot_list[[2]],
#                 plot_list[[3]], plot_list[[4]],
#                 plot_list[[5]], nrow = 1, ncol = 5)
# p2 <- ggarrange(plot_list[[6]], plot_list[[7]],
#                 plot_list[[8]], plot_list[[9]],
#                 plot_list[[10]], nrow = 1, ncol = 5)
# 
# ids <- which(bbs$Time == 1)
# surf <- mba.surf(bbs[ids , c("Longitude","Latitude", "postmedian_z")],
#                  no.X = 200, no.Y = 200, h = 6, m = 1, n = 1,
#                  extend=FALSE)$xyz.est
# surf_rast <- raster(surf)
# us_rast <- raster::extract(surf_rast, usa, cellnumbers = TRUE)
# full_grid <- expand.grid(surf$x, rev(surf$y))
# us_df <- data.frame(full_grid[us_rast[[1]][, "cell"], ],
#                     z = us_rast[[1]][, "value"])
# colnames(us_df) <- c("x", "y", "z")
# us_df <- na.omit(us_df)
# 
# plot_new <- ggplot(usa) + geom_sf() + 
#   coord_sf(xlim = lonlim, ylim = latlim, expand = FALSE) +
#   geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
#   scale_fill_gradientn(colours = col.br(20), breaks = brks, 
#                        label = function(x) sprintf("%.1f", x)) +
#   labs(fill = TeX('$z(s, t)$'), 
#        title = as.character(2009+i)) + 
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.border = element_blank(), 
#         panel.background = element_blank(), 
#         axis.title.x=element_blank(), 
#         axis.title.y=element_blank(), 
#         axis.text.x = element_blank(), 
#         axis.text.y = element_blank(), 
#         axis.ticks = element_blank(), 
#         legend.key.size = unit(0.5, 'cm'),
#         plot.title = element_text(color = "gray40", size = 10,
#                                   hjust = 0.5),
#         aspect.ratio = 0.6)
# legend_new <- get_legend(plot_new)
# 
# finalplot <- ggarrange(p1, p2, nrow = 2, common.legend = TRUE, 
#                        legend = "right", legend.grob = legend_new)
# 
# # system2(command = "pdfcrop",
# #         args    = c("../plots/bbs2010-19_z.pdf",
# #                     "../plots/bbs2010-19_z.pdf")
# # )
# 
