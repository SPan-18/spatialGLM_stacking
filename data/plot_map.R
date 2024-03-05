library("ggplot2")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(ggspatial)
library(basemaps)

world <- ne_countries(scale = "medium", returnclass = "sf")
# class(world)

bbs <- read.csv("BBS21.csv")
bbs <- bbs %>%
  filter(Latitude < 60) %>%
  filter(!(RouteDataID == 6376050)) %>%
  filter(BirdCount < 500)
bbs <- bbs[, c("Latitude", "Longitude", "BirdCount")]

latmin <- min(bbs$Latitude)
latmax <- max(bbs$Latitude)
lonmin <- min(bbs$Longitude)
lonmax <- max(bbs$Longitude)
lat_extra <- 0 #abs(latmax - latmin) * 0.2
lon_extra <- 0 #abs(lonmax - lonmin) * 0.2

latlim <- c(latmin - lat_extra, latmax + lat_extra)
lonlim <- c(lonmin - lon_extra, lonmax + lon_extra)

mybreaks <- quantile(bbs$BirdCount, c(0, 0.2, 0.4, 0.6, 0.8, 1))
blue.red <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
col.br <- colorRampPalette(blue.red)

# ggplot(data = world) +
#   geom_sf(fill = "antiquewhite") +
#   # geom_point(data = bbs, aes(x = Longitude, y = Latitude),
#   #            color = "red") +
#   coord_sf(xlim = lonlim, ylim = latlim, expand = FALSE) +
#   theme_bw() +
#   geom_point(data = bbs, aes(x = Longitude, y = Latitude, 
#                              size = BirdCount, color = BirdCount), 
#              alpha = 1, shape = 20, stroke = FALSE) +
#   # scale_alpha_continuous(trans = "log",
#   #                        range = c(0.1, 0.9), breaks = mybreaks) +
#   scale_color_gradientn(colours = col.br(20), trans = "log",
#                         label = function(x) sprintf("%.0f", x),
#                         name = "Avian\nCount") +
#   scale_size_continuous(trans = "log", 
#                         range = c(1, 5), breaks = mybreaks) +
#   guides(size = 'none') +
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


##using basemaps
# world_2 <- st_transform(world,  crs = st_crs(3857))
# ggplot() +
#   basemap_gglayer(world, map_service = "carto", map_type = "light") +
#   # geom_point(data = bbs, aes(x = Longitude, y = Latitude),
#   #            color = "red") +
#   coord_sf(xlim = lonlim, ylim = latlim, expand = FALSE)
# 
# data(ext_eur)
# # set_defaults(map_service = "osm", map_type = "topographic")
# 
# basemap_magick(ext, map_service = "carto", map_type = "light")
# basemap_magick(world, map_service = "carto", map_type = "light")

### using ggmaps

library(ggmap)
library(viridis)
mapImageData <- get_map(location = c(lon = mean(lonlim), 
                                     lat = mean(latlim)),
                        color = "color", # or bw
                        source = "google",
                        maptype = "satellite",
                        zoom = 3)
                        
ggmap(mapImageData,
      extent = "panel", # "panel" keeps in axes, etc.
      ylab = "Latitude",
      xlab = "Longitude",
      legend = "right",
      darken = 0) +
  geom_point(data = bbs, aes(x = Longitude, y = Latitude, 
                             size = BirdCount, color = BirdCount), 
             alpha = 0.8, shape = 20, stroke = FALSE) +
  # scale_alpha_continuous(trans = "log",
  #                        range = c(0.1, 0.9), breaks = mybreaks) +
  scale_color_gradientn(colours = col.br(20), trans = "log",
                        label = function(x) sprintf("%.0f", x),
                        name = "Raw\nAvian\nCount") +
  scale_size_continuous(trans = "log", 
                        range = c(1.5, 5), breaks = mybreaks) +
  guides(size = 'none') +
  theme_bw() +
  theme(legend.position = c(0.92, 0.18),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size = 8, colour = "antiquewhite"),
        legend.text = element_text(size = 8, colour = "antiquewhite"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        aspect.ratio = 1)


