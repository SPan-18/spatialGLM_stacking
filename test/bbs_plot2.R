rm(list = ls())

bbs <- read.csv("../data/BBS10-19.csv")
bbs <- bbs %>%
  filter(Latitude < 55) %>%
  filter(BirdCount < 500)

post_z <- read.table("bbs_stack/z.txt")

nz <- nrow(bbs)
zcol_ids <- 1*nz + 1:nz
post_zx <- post_z[, zcol_ids]
bbs$plot <- apply(post_zx, 2, median)

library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("sp")
library("raster")
library("ggpubr")

usa <- ne_countries(country = "United States of America",
                    scale = "medium", returnclass = "sf")
# usa_map <- map('usa')
# usa <- st_as_sf(usa_map)
latmin <- min(bbs$Latitude)
latmax <- max(bbs$Latitude)
lonmin <- min(bbs$Longitude)
lonmax <- max(bbs$Longitude)
lat_extra <- abs(latmax - latmin) * 0.00
lon_extra <- abs(lonmax - lonmin) * 0.00
latlim <- c(24.5977, 49.6)
lonlim <- c(-125, -66.5)

blue.red <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
col.br <- colorRampPalette(blue.red)
brks <- c(-0.4, -0.2, 0, 0.2, 0.4)

plot_list_car <- vector(mode = "list", length = 10)

for(i in 1:10){
  ids <- which(bbs$Time == i)
  surf <- mba.surf(bbs[ids , c("Longitude","Latitude", "plot")],
                   no.X = 200, no.Y = 200, h = 6, m = 1, n = 1,
                   extend=FALSE)$xyz.est
  surf_rast <- raster(surf)
  us_rast <- raster::extract(surf_rast, usa, cellnumbers = TRUE)
  full_grid <- expand.grid(surf$x, rev(surf$y))
  us_df <- data.frame(full_grid[us_rast[[1]][, "cell"], ],
                      z = us_rast[[1]][, "value"])
  colnames(us_df) <- c("x", "y", "z")
  us_df <- na.omit(us_df)
  
  plot_list_car[[i]] <- ggplot(usa) +
    geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
    scale_fill_gradientn(colours = col.br(20), breaks = brks, 
                         label = function(x) sprintf("%.1f", x)) +
    geom_sf(data = usa, color = "black", fill = NA, linewidth = 0.05) +
    coord_sf(xlim = lonlim, ylim = latlim, expand = FALSE) +
    labs(fill = TeX('$\\omega_{car}(s, t)$'), 
         title = as.character(2009+i)) + 
    theme_bw() + 
    theme(legend.position = "none", 
          panel.grid = element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.title=element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(color = "gray40", size = 10,
                                    hjust = 0.5, vjust = -1.5),
          plot.margin = unit(c(0,0,0.25,0), "cm"),
          aspect.ratio = 0.75)
}

car1 <- ggarrange(plot_list_car[[1]], plot_list_car[[2]],
                  plot_list_car[[3]], plot_list_car[[4]],
                  plot_list_car[[5]], nrow = 1, ncol = 5)
car2 <- ggarrange(plot_list_car[[6]], plot_list_car[[7]],
                  plot_list_car[[8]], plot_list_car[[9]],
                  plot_list_car[[10]], nrow = 1, ncol = 5)

ids <- which(bbs$Time == 2)
surf <- mba.surf(bbs[ids , c("Longitude","Latitude", "plot")],
                 no.X = 200, no.Y = 200, h = 6, m = 1, n = 1,
                 extend=FALSE)$xyz.est
surf_rast <- raster(surf)
us_rast <- raster::extract(surf_rast, usa, cellnumbers = TRUE)
full_grid <- expand.grid(surf$x, rev(surf$y))
us_df <- data.frame(full_grid[us_rast[[1]][, "cell"], ],
                    z = us_rast[[1]][, "value"])
colnames(us_df) <- c("x", "y", "z")
us_df <- na.omit(us_df)

plot_car_new <- ggplot(usa) + geom_sf() + 
  coord_sf(xlim = lonlim, ylim = latlim, expand = FALSE) +
  geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
  scale_fill_gradientn(colours = col.br(20), breaks = brks, 
                       label = function(x) sprintf("%.1f", x)) +
  labs(fill = TeX('$\\omega_{car}(s, t)$'), 
       title = as.character(2009+i)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_text(color = "gray40", size = 10,
                                  hjust = 0.5),
        aspect.ratio = 0.6)
legend_car <- get_legend(plot_car_new)

carplot <- ggarrange(car1, car2, nrow = 2, common.legend = TRUE, 
                     legend = "right", legend.grob = legend_car)
############################################################################

zcol_ids <- 2*nz + 1:nz
post_zx <- post_z[, zcol_ids]
bbs$plot <- apply(post_zx, 2, median)

plot_list_noise <- vector(mode = "list", length = 10)

for(i in 1:10){
  ids <- which(bbs$Time == i)
  surf <- mba.surf(bbs[ids , c("Longitude","Latitude", "plot")],
                   no.X = 200, no.Y = 200, h = 6, m = 1, n = 1,
                   extend=FALSE)$xyz.est
  surf_rast <- raster(surf)
  us_rast <- raster::extract(surf_rast, usa, cellnumbers = TRUE)
  full_grid <- expand.grid(surf$x, rev(surf$y))
  us_df <- data.frame(full_grid[us_rast[[1]][, "cell"], ],
                      z = us_rast[[1]][, "value"])
  colnames(us_df) <- c("x", "y", "z")
  us_df <- na.omit(us_df)
  
  plot_list_noise[[i]] <- ggplot(usa) +
    geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
    scale_fill_gradientn(colours = col.br(20), breaks = brks, 
                         label = function(x) sprintf("%.1f", x)) +
    geom_sf(data = usa, color = "black", fill = NA, linewidth = 0.05) +
    coord_sf(xlim = lonlim, ylim = latlim, expand = FALSE) +
    labs(fill = TeX('$\\omega_{noise}(s, t)$'), 
         title = as.character(2009+i)) + 
    theme_bw() + 
    theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          axis.title.x=element_blank(), 
          axis.title.y=element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(), 
          plot.title = element_text(color = "gray40", size = 10,
                                    hjust = 0.5, vjust = -1.5),
          plot.margin = unit(c(0,0,0.25,0), "cm"),
          aspect.ratio = 0.75)
}

noise1 <- ggarrange(plot_list_noise[[1]], plot_list_noise[[2]],
                    plot_list_noise[[3]], plot_list_noise[[4]],
                    plot_list_noise[[5]], nrow = 1, ncol = 5)
noise2 <- ggarrange(plot_list_noise[[6]], plot_list_noise[[7]],
                    plot_list_noise[[8]], plot_list_noise[[9]],
                    plot_list_noise[[10]], nrow = 1, ncol = 5)

ids <- which(bbs$Time == 2)
surf <- mba.surf(bbs[ids , c("Longitude","Latitude", "plot")],
                 no.X = 200, no.Y = 200, h = 6, m = 1, n = 1,
                 extend=FALSE)$xyz.est
surf_rast <- raster(surf)
us_rast <- raster::extract(surf_rast, usa, cellnumbers = TRUE)
full_grid <- expand.grid(surf$x, rev(surf$y))
us_df <- data.frame(full_grid[us_rast[[1]][, "cell"], ],
                    z = us_rast[[1]][, "value"])
colnames(us_df) <- c("x", "y", "z")
us_df <- na.omit(us_df)

plot_noise_new <- ggplot(usa) + geom_sf() + 
  coord_sf(xlim = lonlim, ylim = latlim, expand = FALSE) +
  geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
  scale_fill_gradientn(colours = col.br(20), breaks = brks, 
                       label = function(x) sprintf("%.1f", x)) +
  labs(fill = TeX('$\\omega_{noise}(s, t)$'), 
       title = as.character(2009+i)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_text(color = "gray40", size = 18,
                                  hjust = 0.5),
        aspect.ratio = 0.6)
legend_noise <- get_legend(plot_noise_new)

noiseplot <- ggarrange(noise1, noise2, nrow = 2, common.legend = TRUE, 
                       legend = "right", legend.grob = legend_noise)

ggarrange(carplot, noiseplot, nrow = 2, 
          labels = c("(a)", "(b)"), font.label = list(face = "plain"),
          vjust = 1.0, hjust = 0.0)
