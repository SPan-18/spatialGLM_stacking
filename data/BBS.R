library(tidyverse)

YEAR <- 2015

routedat <- read.csv("routes.csv", header = TRUE)
summarydat <- read.csv("MigrantSummary.csv", header = TRUE)
veh <- read.csv("VehicleData.csv", header = TRUE)

dat21 <- summarydat %>%
  filter(Year == YEAR, RPID %in% c(101, 102, 103, 104))

veh21 <- veh %>%
  filter(Year == YEAR) %>%
  rowwise() %>%
  mutate(NCar = sum(c_across(starts_with("Car")), na.rm = T), .keep = "unused") %>%
  mutate(Noise = sum(c_across(starts_with("Noise")), na.rm = T), .keep = "unused") %>%
  dplyr::select(-Year, -RecordedCar)

dat21_merged <- dat21 %>%
  group_by(RouteDataID, CountryNum, StateNum, Route) %>%
  summarize(BirdCount = sum(SpeciesTotal))

dat21_route <- merge(dat21_merged, routedat)

dat21_route <- dat21_route %>%
  dplyr::select(-Active, -Stratum, -BCR, -RouteTypeID, -RouteTypeDetailID)

dat21_route_veh <- merge(dat21_route, veh21)

dat21_final <- dat21_route_veh %>%
  dplyr::select(-RPID) %>%
  dplyr::select(4, 1, 2, 3, 6, 7, 8, 9, 10, 5)

# write.csv(dat21_final, "BBS15.csv")
