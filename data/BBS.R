library(tidyverse)

YEAR <- 2011

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

# Clean duplicate locations

S_dup <- data.frame(table(dat21_final$Latitude))
dup_lats <- S_dup$Var1[S_dup$Freq > 1]
if(length(dup_lats) > 1){
  for(i in 1:length(dup_lats)){
    dup_ids <- which(dat21_final$Latitude == dup_lats[i])
    keep_id <- min(dup_ids)
    dat21_final[keep_id, "NCar"] <- max(dat21_final[dup_ids, "NCar"])
    dat21_final[keep_id, "Noise"] <- max(dat21_final[dup_ids, "Noise"])
    dat21_final[keep_id, "BirdCount"] <- sum(dat21_final[dup_ids, "BirdCount"])
    dat21_final <- dat21_final[- setdiff(dup_ids, keep_id), ]
  }
  cat("Cleaned", length(dup_lats), "duplicate locations..\n")
}

write.csv(dat21_final, paste("bbs/BBS", YEAR-2000, ".csv", sep = ""))
