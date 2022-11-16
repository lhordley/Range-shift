##########################
#### user: Lisbeth Hordley
#### date: July 2021
#### info: Cleaning NMRS data to hectad level for analysis

rm(list = ls())

## load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(sf)
library(ggplot2)
library(plyr)
library(ggpubr)

## read in nmrs data
nmrsdata2 <- fread("Data/NMRS/NMRS_data_range_shift.csv", stringsAsFactors = FALSE) # 670,193 rows

## change column names
colnames(nmrsdata) <- c("Code", "Scientific_name", "Common_name", "Hectad", "Time_period")

## convert 10km gridrefs to lat and lon and Easting and Northing
## after changing elevation to 10x10km this will put everything in 10km hectads for analysis
library(rnrfa)
gridrefs <- as.vector(unique(nmrsdata[,4])) ## 2708 10km gridrefs
lon_lat <- as.data.frame(osg_parse(grid_refs=gridrefs$Hectad, coord_system = "WGS84"))
east_north <- as.data.frame(osg_parse(grid_refs=gridrefs$Hectad))
gridref_latlon_eastnorth <- cbind(gridrefs, lon_lat, east_north)
## convert gridref to east/north - gives easting and northing at 10km scale - used to measure distances between sites later
## merge back into nmrs data
nmrsdata <- merge(nmrsdata, gridref_latlon_eastnorth, by="Hectad")
head(nmrsdata) # 670,193 rows 

## remove Isle of Man records: GB only
# Quick look at data
ggplot(nmrsdata, aes(easting, northing))+
  geom_point()+coord_fixed()

## Need to remove Isle of Man records (GB only) (already removed Northern Ireland by gridref above and no Channel Island records)
# Filter out Isle of Man records ----
sqIM <- nmrsdata[northing < 505000 & northing > 450000 &
                  easting > 200000 & easting < 250000]
ggplot(nmrsdata, aes(easting, northing))+
  geom_point()+coord_fixed()+
  geom_point(aes(easting, northing), data = sqIM, color="red") # looks good

ggplot(nmrsdata[!Hectad %in% sqIM$Hectad], aes(easting, northing))+
  geom_point()+coord_fixed()

## This removes the isle of man hectads
nmrsdata <- nmrsdata[!Hectad %in% sqIM$Hectad] ## 667,521 rows

## save NMRS cleaned data
saveRDS(nmrsdata, "Data/NMRS/NMRS_hectad_cleaned.rds")





















































# 
# ## plot species with migrant populations on a map for Richard to determine which species x hectad combinations to remove
# # Angle-striped Sallow, Great Brocade, Rannoch Looper, Sword Grass, Gold Spangle, Golden-rod Brindle, Red Sword-grass and Scarce Silver Y
# 
# nmrs_migrant <- filter(nmrsdata, (Common_name %in% c("Angle-striped Sallow", "Great Brocade",
#                                                      "Rannoch Looper", "Sword-grass", "Gold Spangle",
#                                                      "Golden-rod Brindle", "Red Sword-grass",
#                                                      "Scarce Silver Y")))
# 
# worldmap = map_data('world')
# migrant_maps <- ggplot() + 
#   geom_polygon(data = worldmap, 
#                aes(x = long, y = lat, group = group), 
#                fill = 'gray90', color = 'black') + 
#   coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#   geom_point(data = nmrs_migrant, 
#              aes(x = as.numeric(lon), 
#                  y = as.numeric(lat)), colour="red", size=0.7) + 
#   theme_void() +
#   theme(title = element_text(size = 12)) +
#   facet_wrap(~Common_name, ncol=4)
# migrant_maps
# ggsave(migrant_maps, file="Maps/Migrant_species_maps.png", height=8, width=10)
# 
# nmrs_migrant <- nmrs_migrant[,c("Common_name", "Hectad", "Time_period")]
# write.csv(nmrs_migrant, file="Data/NMRS/NMRS_migrant_hectads.csv", row.names=FALSE)
# 
# ## re-plot maps with immigrant hectads removed
# nmrs_migrant_exclude <- read.csv("Data/NMRS/NMRS_migrant_hectads_exclude.csv", header=TRUE)
# nmrsdata <- readRDS("Data/NMRS/NMRS_hectad_cleaned.rds")
# 
# nmrs_migrant <- filter(nmrsdata, (Common_name %in% c("Angle-striped Sallow", "Great Brocade",
#                                                      "Rannoch Looper", "Sword-grass", "Gold Spangle",
#                                                      "Golden-rod Brindle", "Red Sword-grass",
#                                                      "Scarce Silver Y")))
# 
# 
# nmrs_migrant <- merge(nmrs_migrant, nmrs_migrant_exclude, by=c("Common_name", "Hectad", "Time_period"))
# ## remove rows where exclude = 1
# nmrs_migrant[is.na(nmrs_migrant)] <- 0
# nmrs_migrant <- nmrs_migrant[!(nmrs_migrant$Exclude==1),]
# 
# nmrs_migrant$Exclude <- as.factor(nmrs_migrant$Exclude)
# worldmap = map_data('world')
# migrant_maps <- ggplot() + 
#   geom_polygon(data = worldmap, 
#                aes(x = long, y = lat, group = group), 
#                fill = 'gray90', color = 'black') + 
#   coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#   geom_point(data = nmrs_migrant, 
#              aes(x = as.numeric(lon), 
#                  y = as.numeric(lat)), colour="indianred2", size=0.7) + 
#   theme_void() +
#   theme(title = element_text(size = 12)) +
#   facet_wrap(~Common_name, ncol=4)
# migrant_maps
# ggsave(migrant_maps, file="Maps/Migrant_species_maps_exclude2.png", height=8, width=10)
# 

# 
# 
# ### explore data to determine which time periods to select for analysis
# ### explore data to determine which time period to use 
# ## plot number of hectads/records per year
# years <- 1970:2016
# worldmap = map_data('world')
# 
# for(i in years) {
#   print(i)
#   # Printing ggplot within for-loop
#   
#   temp_plot <- ggplot() + 
#     geom_polygon(data = worldmap, 
#                  aes(x = long, y = lat, group = group), 
#                  fill = 'gray90', color = 'black') + 
#     coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#     theme_void() + 
#     geom_point(data = nmrsdata[nmrsdata$Year==i,], 
#                aes(x = as.numeric(lon), 
#                    y = as.numeric(lat)), size=1) +
#     ggtitle(i) +
#     theme(title = element_text(size = 12))
#   
#   #ggsave(temp_plot, file=paste0("../../Maps/NMRS_10km_elevation_", i,".png"), width = 20, height = 25, units = "cm")
#   Sys.sleep(2)
# }
# 
# ## massive increase in records from late 90s/early 00s 
# # count the number of recorded hectads and the number of species*hectad records in each year 
# library(plyr)
# nmrsdata$PRESENCE <- 1
# #colnames(nmrsdata)[1] <- "Grid_Square"
# nmrs_rec_hectad_years <- ddply(nmrsdata, .(Hectad, Year, easting, northing), numcolwise(mean)) ## removes species data, just hectads recorded for each year
# 
# recording <- data.frame(Year = numeric(),
#                         Hectads = numeric(),
#                         Records = numeric())
# 
# for (n in years){
#   year_hec <- nmrs_rec_hec_years[which(nmrs_rec_hec_years$Year==n), ] ## NMRS data for year of interest (n)
#   Hectads <- nrow(year_hec) ## number of grid squares recorded in year of interest
#   
#   year_rec <- nmrsdata[which(nmrsdata$Year==n), ] ## find records for unique year sites in nmrs data (i.e. with species data)
#   Records <- nrow(year_rec) ## number of records for year of interest
#   
#   Year <- n
#   
#   out <- data.frame(cbind(Year,Hectads,Records))
#   recording <- rbind(recording,out)
#   
#   
# }
# 
# summary(recording)
# 
# 
# # now plot the number of hectads and records over time
# plot(Grid_Squares ~ Year, data = recording)
# plot(Records ~ Year, data = recording)
# total_grids <- max(recording$Grid_Squares) # 12932 (in 2011)
# 
# recording$Perc <- recording$Grid_Squares*100/total_grids
# 
# plot(Perc ~ Year, data = recording)
# 
# 
# ## need to decide on time periods to investigate
# ## current vs pre-2000
# ## following hickling 2006 - 25 year period between mid-points of time periods
# ## e.g. current is 2012-2016 (mid-point is 2014)
# ## 1989 would be previous mid-point
# ## other time period 1987-1991
# ## longer first time period? By the time we look at heavily recorded hectads, we might have barely any
# ## First TP = 1975-1991 to give more data
# 
# tp1 <- 1975:1991
# tp2 <- 2012:2016
