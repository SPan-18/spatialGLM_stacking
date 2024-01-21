rm(list = ls())

source("../src/runsrc.R")

bbs <- read.csv("../data/BBS21.csv")
bbs <- bbs %>%
  filter(Latitude < 55) %>%
  filter(!(RouteDataID == 6376050)) %>%
  filter(BirdCount < 500)

ids <- 1:nrow(bbs)
y <- as.numeric(bbs[ids, "BirdCount"])
X <- as.matrix(bbs[, c("Longitude", "Latitude", "NCar", "Noise")])
X <- cbind(rep(1, length(y)), X)
S <- as.matrix(bbs[ids, c("Longitude", "Latitude")])

n_postsamp <- 500

mod_list <- create_model_list(G_decay = c(300, 400, 1200), 
                              G_smoothness = c(0.5, 1, 1.5),
                              G_epsilon = c(0.5, 0.75),
                              G_nuxi = 1,
                              G_nubeta = 2.1, G_nuz = 2.1)

m_out <- spGLM_stack(y = y, X = X, S = S, N.samp = n_postsamp,
                     MC.samp = 500,
                     family = "poisson",
                     spCov = "matern",
                     mc.cores = 6,
                     solver = "MOSEK",
                     mod_params_list = mod_list)

postrun_samps <- postrunsampler(m_out, N.samp = n_postsamp)
post_z <- postrun_samps$z
post_beta <- postrun_samps$beta
post_xi <- postrun_samps$xi

print(ci_beta(t(post_beta)))

y.hat <- apply(exp(X %*% post_beta + post_z + post_xi), 2, function(x){rpois(dim(X)[1], x)})
y.hat.mu <- apply(y.hat, 1, median)

bbs$yhat <- log(y.hat.mu)
bbs$postmedian_z <- apply(post_z, 1, median)

library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(sp)
library(raster)
library(viridis)

surf <- mba.surf(bbs[ , c("Longitude","Latitude", "yhat")], 
                 no.X = 200, no.Y = 200, h = 6, m = 1, n = 1, 
                 extend=FALSE)$xyz.est

surf_rast <- raster(surf)
usa <- ne_countries(country = "United States of America", 
                    scale = "medium", returnclass = "sf")
# us_rast <- crop(surf_rast, usa)
us_rast <- raster::extract(surf_rast, usa, cellnumbers = TRUE)
full_grid <- expand.grid(surf$x, rev(surf$y))
us_df <- data.frame(full_grid[us_rast[[1]][, "cell"], ], 
                    z = us_rast[[1]][, "value"])
colnames(us_df) <- c("x", "y", "z")
us_df <- na.omit(us_df)
us_df$z <- exp(us_df$z)

world <- ne_countries(scale = "medium", returnclass = "sf")
latmin <- min(bbs$Latitude)
latmax <- max(bbs$Latitude)
lonmin <- min(bbs$Longitude)
lonmax <- max(bbs$Longitude)
lat_extra <- abs(latmax - latmin) * 0.01
lon_extra <- abs(lonmax - lonmin) * 0.01
latlim <- c(latmin - lat_extra, latmax + lat_extra)
lonlim <- c(lonmin - lon_extra, lonmax + lon_extra)
brks <- quantile(bbs$yhat, c(0, 0.2, 0.4, 0.6, 0.8, 1))
blue.red <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
col.br <- colorRampPalette(blue.red)

ggplot(data = world) +
  geom_sf(fill = "antiquewhite") +
  coord_sf(xlim = lonlim, 
           ylim = latlim, expand = FALSE) +
  # geom_raster(data = surf_df, aes(x = x, y = y, fill = z), alpha = 0.75) +
  geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
  # scale_fill_viridis(option = "plasma", trans = "log", direction = -1, 
  #                    label = function(x) sprintf("%.1f", x)) +
  scale_fill_gradientn(colours = col.br(20), trans = "log",
                       label = function(x) sprintf("%.1f", x)) +
  xlab("Longitude") +
  ylab("Latitude") +
  # labs(fill = latex2exp::TeX('y(s)$')) +
  labs(fill = "Predicted\nAvian\nCount") +
  theme_bw() +
  theme(panel.grid.major = element_line(color = gray(.8), 
                                        linetype = "dotted",
                                        linewidth = 0.25), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = c(0.92, 0.15),
        legend.background = element_blank(),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size = 8),
        legend.title.align = 0.5,
        legend.text = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        aspect.ratio = 1)

# mod_glm <- glm(BirdCount ~ Longitude + Latitude + NCar + Noise, data = bbs,
#                family = poisson(link = "log"))
# summary(mod_glm)



# surf <- mba.surf(bbs[ , c("Longitude","Latitude", "Noise")], 
#                  no.X = 200, no.Y = 200, h = 8, m = 1, n = 1, 
#                  extend=FALSE)$xyz.est
# 
# surf_rast <- raster(surf)
# usa <- ne_countries(country = "United States of America", 
#                     scale = "medium", returnclass = "sf")
# # us_rast <- crop(surf_rast, usa)
# us_rast <- raster::extract(surf_rast, usa, cellnumbers = TRUE)
# full_grid <- expand.grid(surf$x, rev(surf$y))
# us_df <- data.frame(full_grid[us_rast[[1]][, "cell"], ], 
#                     z = us_rast[[1]][, "value"])
# colnames(us_df) <- c("x", "y", "z")
# us_df <- na.omit(us_df)
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# brks <- quantile(bbs$yhat, c(0, 0.2, 0.4, 0.6, 0.8, 1))
# 
# blue.red <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# col.br <- colorRampPalette(blue.red)
# 
# ggplot(data = world) +
#   geom_sf(fill = "white") +
#   coord_sf(xlim = c(-149.40764, -39.89025), 
#            ylim = c(17.06898, 76.91902), expand = FALSE) +
#   # geom_raster(data = surf_df, aes(x = x, y = y, fill = z), alpha = 0.75) +
#   geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
#   # scale_fill_viridis(option = "plasma", trans = "log", direction = -1, 
#   #                    label = function(x) sprintf("%.1f", x)) +
#   scale_fill_gradientn(colours = col.br(20),
#                        label = function(x) sprintf("%.1f", x)) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   # labs(fill = latex2exp::TeX('y(s)$')) +
#   labs(fill = "Predicted\nAvian\nCount") +
#   theme_bw() +
#   theme(panel.grid.major = element_line(color = gray(.8), 
#                                         linetype = "dotted",
#                                         linewidth = 0.25), 
#         panel.background = element_rect(fill = "aliceblue"),
#         legend.position = c(0.92, 0.15),
#         legend.background = element_blank(),
#         legend.key.size = unit(0.3, 'cm'),
#         legend.title = element_text(size = 8),
#         legend.title.align = 0.5,
#         legend.text = element_text(size = 8),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         aspect.ratio = 1)




