
routedat <- read.csv("routes.csv", header = TRUE)
summarydat <- read.csv("MigrantSummary.csv", header = TRUE)

dat21 <- summarydat %>% 
  filter(Year == 2021)

# summary_CA2021 <- subset(summarydat, Year == '2021'& CountryNum == '124')
#                          # select = c('Route', 'SpeciesTotal'))
# summary_CA2021$lat <- NA
# summary_CA2021$lon <- NA
# 
# route_CA <- subset(routedat, CountryNum == '124')
#                    # select = c('Route', 'Latitude', 'Longitude'))
# 
# for(i in 1:nrow(summary_CA2021)){
#   rowid <- which(route_CA$Route == summary_CA2021$Route[i])
#   summary_CA2021$lat[i] <- route_CA[rowid, "Latitude"]
#   summary_CA2021$lon[i] <- route_CA[rowid, "Longitude"]
# }
# 
# summary_CA2021 <- subset(summary_CA2021, 
#                          select = c("lat", "long", 'SpeciesTotal'))


library(leaflet)
library(sp)
# route_coords <- routedat[, c("Latitude", "Longitude")]
# colnames(route_coords) <- c("lat", "lon")
# coords <- route_coords
# temp <- rnorm(nrow(route_coords))

# coordinates(coords) <- ~lon+lat
# proj4string(coords) <- "+proj=longlat +datum=WGS84"
# coords.utm <- spTransform(coords, CRS("+proj=utm +zone=13 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
# coords.utm <- coordinates(coords.utm)/1000

base.map <- leaflet(width="100%") %>%
  addProviderTiles("Esri.WorldImagery", group="Satellite") %>%
  addLayersControl(
    baseGroup = c("Satellite"),
    options = layersControlOptions(collapsed = FALSE)
  )

base.map %>%
  addCircleMarkers(lng = dat21_route_veh[,"Longitude"], 
                   lat = dat21_route_veh[,"Latitude"], 
                   col = "red", stroke = FALSE, 
                   radius = 3, fillOpacity = 0.9)
