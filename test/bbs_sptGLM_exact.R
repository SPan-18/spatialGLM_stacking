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

n_postsamp <- 1000

# mod_list <- create_model_list(G_decay = c(400, 800), 
#                               G_smoothness = c(0.5, 1),
#                               G_timedecay = c(0.9),
#                               G_epsilon = c(0.5),
#                               G_nuxi = 1,
#                               G_nubeta = 2.1, G_nuz = 2.1)

params <- list(phi = 800, phi_t = 0.9, nu_matern = 0.5,
               nu_xi = 1, nu_beta = 2.1, nu_z = 2.1,
               alpha_epsilon = 0.5)

set.seed(1729)
t1 <- Sys.time()
mod_out <- sptGLM_exact(y = y, X = X, S = S, time = time,
                        N.samp = n_postsamp,
                        family = "poisson",
                        mod_params = params)
t2 <- Sys.time()
cat("Time took:", t2-t1, units(t2-t1), ".\n")

post_beta <- mod_out$beta
post_z <- mod_out$z
post_xi <- mod_out$xi

write.table(post_beta,
            file = "bbs_stack/beta.txt",
            col.names = FALSE, row.names = FALSE)
write.table(post_z,
            file = "bbs_stack/z.txt",
            col.names = FALSE, row.names = FALSE)
write.table(post_xi,
            file = "bbs_stack/xi.txt",
            col.names = FALSE, row.names = FALSE)


# print(ci_beta(post_beta))

# y.hat <- t(apply(exp(tcrossprod(X, post_beta) +
#                      tcrossprod(makeG(X), post_z) + t(post_xi)),
#                 1, function(x){rpois(n_postsamp, x)}))
# y.hat.mu <- apply(y.hat, 1, function(x) quantile(x, 0.5))
# 
# bbs$yhat <- log(y.hat.mu)

# nz <- length(y)
# zcol_ids <- 2*nz + 1:nz
# post_zx <- mod_out$z[, zcol_ids]
# bbs$yhat <- apply(post_zx, 2, median)
# #
# #
# library("sf")
# library("rnaturalearth")
# library("rnaturalearthdata")
# library(sp)
# library(raster)
# library(viridis)
# library(mapdata)
# 
# plot_list <- vector(mode = "list", length = 10)
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# latmin <- min(bbs$Latitude)
# latmax <- max(bbs$Latitude)
# lonmin <- min(bbs$Longitude)
# lonmax <- max(bbs$Longitude)
# lat_extra <- abs(latmax - latmin) * 0.01
# lon_extra <- abs(lonmax - lonmin) * 0.01
# latlim <- c(latmin - lat_extra, latmax + lat_extra)
# lonlim <- c(lonmin - lon_extra, lonmax + lon_extra)
# blue.red <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# col.br <- colorRampPalette(blue.red)
# brks <- c(-0.4, -0.2, 0, 0.2, 0.4)
# 
# for(i in 1:10){
#   ids <- which(bbs$Time == i)
#   surf <- mba.surf(bbs[ids , c("Longitude","Latitude", "yhat")],
#                    no.X = 200, no.Y = 200, h = 6, m = 1, n = 1,
#                    extend=FALSE)$xyz.est
#   surf_rast <- raster(surf)
#   usa <- ne_countries(country = "United States of America",
#                       scale = "medium", returnclass = "sf")
#   # us_rast <- crop(surf_rast, usa)
#   us_rast <- raster::extract(surf_rast, usa, cellnumbers = TRUE)
#   full_grid <- expand.grid(surf$x, rev(surf$y))
#   us_df <- data.frame(full_grid[us_rast[[1]][, "cell"], ],
#                       z = us_rast[[1]][, "value"])
#   colnames(us_df) <- c("x", "y", "z")
#   us_df <- na.omit(us_df)
#   us_df$z <- exp(us_df$z)
#   # brks <- quantile(us_df$z, c(0, 0.4, 0.8, 0.95, 0.995, 1))
#   # brks <- quantile(us_df$z, c(0, 0.1, 0.9, 1))
# 
#   # plot_list[[i]] <- ggplot(data = world) +
#   #   geom_sf(fill = "antiquewhite") +
#   #   coord_sf(xlim = lonlim,
#   #            ylim = c(25, 49.5), expand = FALSE) +
#   #   geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
#   #   # scale_fill_viridis(option = "plasma", trans = "log", direction = -1,
#   #   #                    label = function(x) sprintf("%.1f", x)) +
#   #   scale_fill_gradientn(colours = col.br(20), #trans = "log",
#   #                        breaks = brks,
#   #                        label = function(x) sprintf("%.1f", x)) +
#   #   xlab("Longitude") +
#   #   ylab("Latitude") +
#   #   # labs(fill = latex2exp::TeX('y(s)$')) +
#   #   # labs(fill = "Predicted\nCount",
#   #   #      title = as.character(2009+i)) +
#   #   labs(fill = TeX('$\\omega_{noise}(s, t)$'),
#   #        title = as.character(2009+i)) +
#   #   theme_bw() +
#   #   theme(panel.grid.major = element_line(color = gray(.8),
#   #                                         linetype = "dotted",
#   #                                         linewidth = 0.25),
#   #         panel.background = element_rect(fill = "aliceblue"),
#   #         legend.position = c(0.92, 0.15),
#   #         legend.background = element_blank(),
#   #         legend.key.size = unit(1, 'cm'),
#   #         legend.title = element_text(size = 14),
#   #         legend.title.align = 0.5,
#   #         legend.text = element_text(size = 12),
#   #         axis.title.x=element_blank(),
#   #         axis.title.y=element_blank(),
#   #         plot.title = element_text(color = "gray40", size = 18,
#   #                                   hjust = 0.5),
#   #         aspect.ratio = 0.8)
# 
#   plot_list[[i]] <- ggplot(usa) +
#     geom_sf() +
#     coord_sf(xlim = lonlim,
#              ylim = c(25, 49.5), expand = FALSE) +
#     geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
#     scale_fill_gradientn(colours = col.br(20), #trans = "log",
#                          breaks = brks,
#                          label = function(x) sprintf("%.1f", x)) +
#     labs(fill = TeX('$\\omega_{noise}(s, t)$'),
#          title = as.character(2009+i)) +
#     theme_bw() +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           panel.background = element_blank(),
#           legend.position = c(0.92, 0.15),
#           legend.background = element_blank(),
#           legend.key.size = unit(1, 'cm'),
#           legend.title = element_text(size = 14),
#           legend.title.align = 0.5,
#           legend.text = element_text(size = 12),
#           axis.title.x=element_blank(),
#           axis.title.y=element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank(),
#           axis.ticks = element_blank(),
#           plot.title = element_text(color = "gray40", size = 18,
#                                     hjust = 0.5),
#           aspect.ratio = 0.67)
# }
# 
# # gridExtra::grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]],
# #                         plot_list[[4]], plot_list[[5]], plot_list[[6]],
# #                         plot_list[[7]], plot_list[[8]], plot_list[[9]],
# #                         ncol = 3)
# 
# library(ggpubr)
# ggarrange(plot_list[[2]], plot_list[[3]],
#           plot_list[[4]], plot_list[[5]], plot_list[[6]],
#           plot_list[[7]], plot_list[[8]], plot_list[[9]],
#           plot_list[[10]],
#           ncol = 3, nrow = 3, common.legend = TRUE, legend="bottom")
# 
