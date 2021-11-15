##########################
#### user: Lisbeth Hordley
#### date: October 2021
#### info: Add temperature data for NMRS sites

rm(list = ls())

## load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(sf)
library(ggplot2)
library(plyr)
library(here)

##########################################################################################################
########### EARLY TIME PERIOD: 1975 - 1991 USED TO SELECT UPLAND/NORTHERN SPCEIES FOR ANALYSIS ###########
##########################################################################################################

## Extract annual temperature data & merge with NMRS data

## read in temperature data (5km scale from metoffice CP18, annual temperature for years 1975:1991)
## https://catalogue.ceda.ac.uk/uuid/f2da35c56afb4fa6aebf44094b65dff3 
library(ncdf4)
library(raster)
files <- list.files(path="Data/UKCP_tas_annual/Early_TP_75_91/",
                    pattern='*.nc',full.names=TRUE)

s <- raster::stack(files)
s[s %in% 9.969209968386869047442886268468442020e+36] <- NA ## change these values to NA

s2 = raster(extent(s), resolution = 10000, crs = crs(s)) ## change resolution to 10x10km
s3 = resample(s, s2)

mean <- stackApply(s3, indices =  rep(1,nlayers(s3)), fun = "mean", na.rm = T) ## take mean across the 17 years/files
plot(mean) # looks good

mean_newproj <- projectRaster(mean,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ## reproject to correct projection
plot(mean_newproj, axes = T) ## mean temperature averaged across years 1975-1991 in 10x10km squares 

##

mean_temp <- as.data.frame(rasterToPoints(mean_newproj)) ## get temp values as dataframe (2929 rows)
colnames(mean_temp) <- c("lon", "lat", "temperature")
## save temp data
write.csv(mean_temp, file="Data/UKCP_tas_annual/Early_TP_75_91/mean_temp_75_91.csv", row.names=FALSE)
mean_temp <- read.csv("Data/UKCP_tas_annual/Early_TP_75_91/mean_temp_75_91.csv", header=TRUE)

## better plot of temperature data
library(ggplot2)
worldmap = map_data('world')
temp_map <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = nmrsdata, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), size=1) + 
  theme_void() +
  theme(title = element_text(size = 12))
temp_map

## now match these with NMRS data
nmrsdata <- readRDS("Data/NMRS/NMRS_hectad_elevation.rds")
## add 5km to easting and northing values to move record to centre of hectad
## this will ensure we are taking the temperature which represents middle of hectad
nmrsdata$easting_centre <- nmrsdata$easting + 5000
nmrsdata$northing_centre <- nmrsdata$northing + 5000

## convert easting + northing to lat + lon
wgs84 = "+init=epsg:4326"
bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 
+ellps=airy +datum=OSGB36 +units=m +no_defs'

ConvertCoordinates <- function(easting,northing) {
  out = cbind(easting,northing)
  mask = !is.na(easting)
  sp <-  sp::spTransform(sp::SpatialPoints(list(easting[mask],northing[mask]),proj4string=sp::CRS(bng)),sp::CRS(wgs84))
  out[mask,]=sp@coords
  out
}

x <- ConvertCoordinates(nmrsdata$easting_centre, nmrsdata$northing_centre)
colnames(x) <- c("lon_centre","lat_centre")
nmrsdata <- cbind(nmrsdata, x) ## nmrsdata now has lat and lon values for centre of hectad to match to nearest temperature hectads

## now match centre lat/lon NMRS values with nearest mean_temp lat/lon values
nmrs_latlon <- nmrsdata[,c(13:14)]
nmrs_latlon <- data.frame(unique(nmrs_latlon)) ## 2719 rows

## need to source Find Closest Rcpp function first!
Rcpp::sourceCpp("Scripts/Functions/Find_Closest_Rcpp_Function.cpp")

find_closest <- function(lat1, long1, lat2, long2) {
  
  toRad <- pi / 180
  lat1  <- lat1  * toRad
  long1 <- long1 * toRad
  lat2  <- lat2  * toRad
  long2 <- long2 * toRad
  
  ord1  <- order(lat1)
  rank1 <- match(seq_along(lat1), ord1)
  ord2  <- order(lat2)
  
  ind <- find_closest_point(lat1[ord1], long1[ord1], lat2[ord2], long2[ord2])
  
  ord2[ind + 1][rank1]
}


result <- find_closest(nmrs_latlon[, 2], nmrs_latlon[, 1], 
                       mean_temp[, 2], mean_temp[, 1])

result ### works super fast!!! assume numbers are the corresponding row numbers that each lat_lon value matches in temp_latlon?
## 62539 results, but 59967 unique values - some temperature sites have matched with multiple NMRS sites (mainly on the coast & Spurn Head)
## add row numbers to temp_latlon
nmrs_latlon$row_numbers <- rownames(nmrs_latlon)
nmrs_latlon$temp_row <- result
mean_temp$row_numbers <- rownames(mean_temp)

lat_lon_final <- merge(nmrs_latlon, mean_temp, by.x="temp_row", by.y="row_numbers")
lat_lon_final <- lat_lon_final[,-c(1,4:6)]
colnames(lat_lon_final) <- c("lon_centre", "lat_centre", "temperature")

## map of mean air temp at NMRS 10km sqaures/points
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = lat_lon_final, 
             aes(x = as.numeric(lon_centre), 
                 y = as.numeric(lat_centre), colour = temperature), size=1) + 
  scale_color_viridis_c(name="Mean temperature 1975-1991") + 
  theme_void() +
  theme(title = element_text(size = 12))

## merge in with NMRS data
nmrsdata_temp <- merge(nmrsdata, lat_lon_final, by=c("lon_centre", "lat_centre"))
## need to subset by year (temperature values only refer to 1975-1991 so occupied hectads need to reflect this too)
## then remove year column and unique - so we only have a list of unique hectads for each species
nmrsdata_temp_early <- nmrsdata_temp[nmrsdata_temp$Year>=1975 & nmrsdata_temp$Year<=1991, ]
unique(nmrsdata_temp_early$Year) ## 17 years
## now remove year and unique dataset - this means that when a species is recorded in the same hectad
## in multiple years, it's not counted in our calculation of the median
nmrsdata_temp_early$Year <- NULL
nmrsdata_temp_early <- unique(nmrsdata_temp_early) ## 165005 rows
## save data: NMRS data for all species & all hectads with mean annual temperature 
saveRDS(nmrsdata_temp_early, file="Data/NMRS/NMRS_hectad_temperature_TP1.rds") ## 165005 rows

########################################################################################################
## PLOT SPECIES TEMPERATURE RANGES (MEDIAN & 75TH PERCENTILE)
## These are given to Richard to determine which cool-adapted species to select
rm(list = ls())

nmrsdata_temp_early <- readRDS("../../Range shift/Data/NMRS/NMRS_hectad_temperature_TP1.rds") ## NMRS data with temperature 
migrant_hectads_early <- read.csv("../../Range shift/Data/NMRS/migrant_hectad_exclusion_TP1.csv", header=TRUE) ## species x hectad combinations to remove due to immigrant populations
## provided by Richard in Teams 20/10/21

# ## save hectad x species data for species with migrant populations (Richard used this to create migrant_hectads_exclusion_early)
# ## subset by species of interest
# nmrsdata_temp_early_filter <- filter(nmrsdata_temp_early, (Common_Name=="Sword-grass") | (Common_Name=="Great Brocade") | 
#                                        (Common_Name=="Angle-striped Sallow"))
# nmrsdata_temp_late_filter <- filter(nmrsdata_temp_late, (Common_Name=="Sword-grass") | (Common_Name=="Great Brocade") | 
#                                       (Common_Name=="Angle-striped Sallow") | Common_Name=="Rannoch Looper")
# ## only take common name and hectad
# nmrsdata_temp_early_filter <- nmrsdata_temp_early_filter[,c(5,3)]
# nmrsdata_temp_late_filter <- nmrsdata_temp_late_filter[,c(5,3)]
# ## save files
# write.csv(nmrsdata_temp_early_filter, file="../../Range shift/Data/NMRS/migrant_hectad_TP1.csv", row.names=FALSE)
# write.csv(nmrsdata_temp_late_filter, file="../../Range shift/Data/NMRS/migrant_hectad_TP2.csv", row.names=FALSE)

## merge migrant exclusion data with nmrs_temp data
## remove hectad/species combinations where exclude == 1
migrant_hectads_early[is.na(migrant_hectads_early)] <- 0
nmrsdata_temp_early <- merge(nmrsdata_temp_early, migrant_hectads_early, by=c("Common_Name", "Hectad"), all=TRUE)
nmrsdata_temp_early[is.na(nmrsdata_temp_early)] <- 0
nmrsdata_temp_early <- nmrsdata_temp_early[nmrsdata_temp_early$Exclude !=1, ] ## removes 84 rows (now 165390)

detach(package:plyr)
library(dplyr)
library(ggplot2)
summary_nmrs_temp <- nmrsdata_temp_early %>% 
  group_by(Common_Name) %>% 
  summarise(lower = min(temperature), upper = max(temperature), p = median(temperature)) ## 796 species
summary_nmrs_temp$lower <- as.numeric(summary_nmrs_temp$lower)
summary_nmrs_temp$upper <- as.numeric(summary_nmrs_temp$upper)
summary_nmrs_temp$p <- as.numeric(summary_nmrs_temp$p)

summary_nmrs_temp <- summary_nmrs_temp[order(summary_nmrs_temp$p, decreasing = TRUE),]  
library(Hmisc)
summary_nmrs_temp$groups<-as.numeric(cut2(summary_nmrs_temp$p, g=12))
## these groups are used to put groups into nmrs_temp to plot them more easily
sp_groups <- summary_nmrs_temp[,c(1,5)]

library(Hmisc)
## use groups from above to put each species into groups to make plots
nmrsdata_temp_early <- merge(nmrsdata_temp_early, sp_groups, by="Common_Name", all=TRUE)
group <- unique(nmrsdata_temp_early$groups) ## 12 groups

for(i in group) {
  print(i)
  # Printing ggplot within for-loop
  
  # temp_plot <- ggplot(data = summary_nmrs_temp[summary_nmrs_temp$groups==i,], mapping = aes(x = reorder(Common_Name, -p), y = p)) +
  #   geom_pointrange(mapping = aes(ymin = lower, ymax = upper)) +
  #   coord_flip() +
  #   labs(x="Common name", y="Median temperature") +
  #   ylim(3,12) +
  #   theme_classic()
  temp_plot <- ggplot(data = nmrsdata_temp_early[nmrsdata_temp_early$groups==i,], mapping = aes(x = reorder(Common_Name,-temperature, FUN=median), y = temperature)) +
    geom_boxplot() +
    #geom_hline(yintercept=9, linetype="dotted") +
    coord_flip() +
    labs(x="Common name", y="Median temperature") +
    scale_y_continuous(breaks = seq(3, 11, by = 1)) +
    theme_light()
  
  ggsave(temp_plot, file=paste0("../../Range shift/Data/UKCP_tas_annual/Early_TP_75_91/NMRS_temp_range_box_", i,".png"), width = 20, height = 25, units = "cm")
  Sys.sleep(2)
}

## calculate median and 75th percentile for each species for Richard
temp_ranges <- nmrsdata_temp_early %>%  group_by(Common_Name) %>% 
  summarise(median_temperature = median(temperature), percentile_75th = quantile(temperature, probs = 0.75))
## order by median ascending
temp_ranges <- arrange(temp_ranges, median_temperature)
## save file
write.csv(temp_ranges, file="../../Range shift/Data/NMRS_temp_median_percentile.csv", row.names=FALSE)
## Boxplots and percentile list given to Richard to select upland species

###################################################################################################

## create upland/northern species list using Richard's classification
## Provided on Teams on 20/10/21

upland_species <- read.csv("../../Range shift/Data/NMRS/NMRS_temp_median_percentile_selection.csv", header=TRUE)
upland_species[is.na(upland_species)] <- 0 # replace NAs with 0 (0 = include the species as northern/upland)
upland_species_final <- upland_species[upland_species$exclude !=1, ] ## 73 species
upland_species_final <- upland_species_final[,1]
upland <- upland_species_final[,1]
nmrsdata_temp_early_upland <- nmrsdata_temp_early[which(nmrsdata_temp_early$Common_Name %in% upland), ] 
## 5390 unique upland species x hectads between 1975-1991

## plot map - make sure there are no hectads in lowland areas
## map of mean air temp at NMRS 10km sqaures/points
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = nmrsdata_temp_early_upland, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), size=1) + 
  theme_void() +
  theme(title = element_text(size = 12))







##########################################################################################################
######### LATE TIME PERIOD: 2012 - 2016 COMPARE MEDIAN/75TH PERCENTILE TEMPERATURE TO EARLY TP ###########
##########################################################################################################

## Extract annual temperature data & merge with NMRS data

## read in temperature data (5km scale from metoffice CP18, annual temperature for years 2012:2016)
## https://catalogue.ceda.ac.uk/uuid/f2da35c56afb4fa6aebf44094b65dff3 
library(ncdf4)
library(raster)
files <- list.files(path="Data/UKCP_tas_annual/Late_TP_12_16/",
                    pattern='*.nc',full.names=TRUE)

s <- raster::stack(files)
s[s %in% 9.969209968386869047442886268468442020e+36] <- NA ## change these values to NA

s2 = raster(extent(s), resolution = 10000, crs = crs(s)) ## change resolution to 10x10km
s3 = resample(s, s2)

mean <- stackApply(s3, indices =  rep(1,nlayers(s3)), fun = "mean", na.rm = T) ## take mean across the 5 years/files
plot(mean) # looks good

mean_newproj <- projectRaster(mean,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ## reproject to correct projection
plot(mean_newproj, axes = T) ## mean temperature averaged across years 1975-1991 in 10x10km squares 

##

mean_temp <- as.data.frame(rasterToPoints(mean_newproj)) ## get temp values as dataframe (2929 rows)
colnames(mean_temp) <- c("lon", "lat", "temperature")
## save file
write.csv(mean_temp, file="Data/UKCP_tas_annual/Late_TP_12_16/mean_temp_12_16.csv", row.names=FALSE)
mean_temp <- read.csv("Data/UKCP_tas_annual/Late_TP_12_16/mean_temp_12_16.csv", header=TRUE)

## better plot of temperature data
library(ggplot2)
worldmap = map_data('world')
temp_map <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = mean_temp, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour = temperature), size=1) + 
  scale_color_viridis_c(name="Mean air temperature") + 
  theme_void() +
  theme(title = element_text(size = 12))
temp_map

## now match these with NMRS data
nmrsdata <- readRDS("Data/NMRS/NMRS_hectad_elevation.rds")
## add 5km to easting and northing values to move record to centre of hectad
## this will ensure we are taking the temperature which represents middle of hectad
nmrsdata$easting_centre <- nmrsdata$easting + 5000
nmrsdata$northing_centre <- nmrsdata$northing + 5000

## convert easting + northing to lat + lon
wgs84 = "+init=epsg:4326"
bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 
+ellps=airy +datum=OSGB36 +units=m +no_defs'

ConvertCoordinates <- function(easting,northing) {
  out = cbind(easting,northing)
  mask = !is.na(easting)
  sp <-  sp::spTransform(sp::SpatialPoints(list(easting[mask],northing[mask]),proj4string=sp::CRS(bng)),sp::CRS(wgs84))
  out[mask,]=sp@coords
  out
}

x <- ConvertCoordinates(nmrsdata$easting_centre, nmrsdata$northing_centre)
colnames(x) <- c("lon_centre","lat_centre")
nmrsdata <- cbind(nmrsdata, x) ## nmrsdata now has lat and lon values for centre of hectad to match to nearest temperature hectads

## now match centre lat/lon NMRS values with nearest mean_temp lat/lon values
nmrs_latlon <- nmrsdata[,c(13:14)]
nmrs_latlon <- data.frame(unique(nmrs_latlon)) ## 2719 rows

## need to source Find Closest Rcpp function first!
Rcpp::sourceCpp("Scripts/Functions/Find_Closest_Rcpp_Function.cpp")

find_closest <- function(lat1, long1, lat2, long2) {
  
  toRad <- pi / 180
  lat1  <- lat1  * toRad
  long1 <- long1 * toRad
  lat2  <- lat2  * toRad
  long2 <- long2 * toRad
  
  ord1  <- order(lat1)
  rank1 <- match(seq_along(lat1), ord1)
  ord2  <- order(lat2)
  
  ind <- find_closest_point(lat1[ord1], long1[ord1], lat2[ord2], long2[ord2])
  
  ord2[ind + 1][rank1]
}


result <- find_closest(nmrs_latlon[, 2], nmrs_latlon[, 1], 
                       mean_temp[, 2], mean_temp[, 1])

result ### works super fast!!! assume numbers are the corresponding row numbers that each lat_lon value matches in temp_latlon?
## 62539 results, but 59967 unique values - some temperature sites have matched with multiple NMRS sites (mainly on the coast & Spurn Head)
## add row numbers to temp_latlon
nmrs_latlon$row_numbers <- rownames(nmrs_latlon)
nmrs_latlon$temp_row <- result
mean_temp$row_numbers <- rownames(mean_temp)

lat_lon_final <- merge(nmrs_latlon, mean_temp, by.x="temp_row", by.y="row_numbers")
lat_lon_final <- lat_lon_final[,-c(1,4:6)]
colnames(lat_lon_final) <- c("lon_centre", "lat_centre", "temperature")

## map of mean air temp at NMRS 10km sqaures/points
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = lat_lon_final, 
             aes(x = as.numeric(lon_centre), 
                 y = as.numeric(lat_centre), colour = temperature), size=1) + 
  scale_color_viridis_c(name="Mean temperature 1975-1991") + 
  theme_void() +
  theme(title = element_text(size = 12))

## merge in with NMRS data
nmrsdata_temp <- merge(nmrsdata, lat_lon_final, by=c("lon_centre", "lat_centre"))
## need to subset by year (temperature values only refer to 1975-1991 so occupied hectads need to reflect this too)
## then remove year column and unique - so we only have a list of unique hectads for each species
nmrsdata_temp_late <- nmrsdata_temp[nmrsdata_temp$Year>=2012 & nmrsdata_temp$Year<=2016, ]
unique(nmrsdata_temp_late$Year) ## 5 years
## now remove year and unique dataset - this means that when a species is recorded in the same hectad
## in multiple years, it's not counted in our calculation of the median
nmrsdata_temp_late$Year <- NULL
nmrsdata_temp_late <- unique(nmrsdata_temp_late) ## 388890 rows
## save data: NMRS data for all species & unique hectads with mean annual temperature for years 2012-2016
saveRDS(nmrsdata_temp_late, file="Data/NMRS/NMRS_hectad_temperature_TP2.rds")






















































# 
# ################
# ## 10km scale ##
# ################
# 
# ## Summer temperature
# 
# ## read in temperature data (5km scale from metoffice CP18, monthly temperature for years 1975:1991)
# library(ncdf4)
# library(raster)
# files <- list.files(path="../../Data/UKCP_tas_monthly/",
#                     pattern='*.nc',full.names=TRUE)
# 
# s <- raster::stack(files)
# #s[s %in% 9.969209968386869047442886268468442020e+36] <- NA ## change these values to NA
# 
# 
# s2 = raster(extent(s), resolution = 10000, crs = crs(s)) ## change resolution to 10x10km
# s3 = resample(s, s2)
# 
# mean_newproj <- projectRaster(s3,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ## reproject to correct projection
# plot(mean_newproj, axes = T) ## monthly temp for each of the 17 years (haven't taken mean yet)
# 
# ## June, July & August = Summer months
# mean_temp <- as.data.frame(rasterToPoints(mean_newproj)) ## get temp values as dataframe (2929 rows, 206 columns (2 of them are lat/lon))
# ## need to rep years (1975:1991) and rep months (x12) to make this into a long dataframe
# library(reshape2)
# mean_temp_final <- melt(data=mean_temp, id.vars = c("x", "y"), variable.name = "dataset", value.name = "mean_temp")
# months = rep(c("Jan","Feb","March","April","May","Jun","Jul","Aug","Sept","Oct",
#                "Nov", "Dec"), times=17, each=2929)
# years <- rep(1975:1991, each=35148)
# mean_temp_final$month <- months
# mean_temp_final$year <- years
# ## remove dataset
# mean_temp_final$dataset <- NULL
# colnames(mean_temp_final)[1:2] <- c("lon","lat")
# 
# ## now only keep summer months
# mean_temp_final <- subset(mean_temp_final, (month %in% c("Jun","Jul","Aug")))
# ## now take mean across summer months and across years (left with one temperature value for each hectad)
# library(dplyr)
# mean_temp_final2 <- mean_temp_final %>% group_by(lon, lat) %>% dplyr::summarize(summer_mean_temp=mean(mean_temp))
# mean_temp_final2 <- as.data.frame(mean_temp_final2)
# ## better plot of temperature data
# library(ggplot2)
# worldmap = map_data('world')
# summer_temp_map <- ggplot() + 
#   geom_polygon(data = worldmap, 
#                aes(x = long, y = lat, group = group), 
#                fill = 'gray90', color = 'black') + 
#   coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#   geom_point(data = mean_temp_final2, 
#              aes(x = as.numeric(lon), 
#                  y = as.numeric(lat), colour = summer_mean_temp), size=1) + 
#   scale_color_viridis_c(name="Mean summer temperature") + 
#   theme_void() +
#   theme(title = element_text(size = 12))
# summer_temp_map
# 
# ## now match these with NMRS data
# nmrsdata <- readRDS("../../Data/NMRS/NMRS_hectad_cleaned.rds")
# ## add 5km to easting and northing values to move record to centre of hectad
# ## this will ensure we are taking the temperature which represents middle of hectad
# nmrsdata$easting_centre <- nmrsdata$easting + 5000
# nmrsdata$northing_centre <- nmrsdata$northing + 5000
# 
# wgs84 = "+init=epsg:4326"
# bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 
# +ellps=airy +datum=OSGB36 +units=m +no_defs'
# 
# ConvertCoordinates <- function(easting,northing) {
#   out = cbind(easting,northing)
#   mask = !is.na(easting)
#   sp <-  sp::spTransform(sp::SpatialPoints(list(easting[mask],northing[mask]),proj4string=sp::CRS(bng)),sp::CRS(wgs84))
#   out[mask,]=sp@coords
#   out
# }
# 
# x <- ConvertCoordinates(nmrsdata$easting_centre, nmrsdata$northing_centre)
# colnames(x) <- c("lon_centre","lat_centre")
# nmrsdata <- cbind(nmrsdata, x)
# 
# ## now match centre lat/lon NMRS values with nearest mean_temp lat/lon values
# nmrs_latlon <- nmrsdata[,c(11:12)]
# nmrs_latlon <- data.frame(unique(nmrs_latlon))
# 
# ## need to source Find Closest Rcpp function first!
# 
# find_closest <- function(lat1, long1, lat2, long2) {
#   
#   toRad <- pi / 180
#   lat1  <- lat1  * toRad
#   long1 <- long1 * toRad
#   lat2  <- lat2  * toRad
#   long2 <- long2 * toRad
#   
#   ord1  <- order(lat1)
#   rank1 <- match(seq_along(lat1), ord1)
#   ord2  <- order(lat2)
#   
#   ind <- find_closest_point(lat1[ord1], long1[ord1], lat2[ord2], long2[ord2])
#   
#   ord2[ind + 1][rank1]
# }
# 
# 
# result <- find_closest(nmrs_latlon[, 2], nmrs_latlon[, 1], 
#                        mean_temp_final2[, 2], mean_temp_final2[, 1])
# 
# result ### works super fast!!! assume numbers are the corresponding row numbers that each lat_lon value matches in temp_latlon?
# ## add row numbers to temp_latlon
# nmrs_latlon$row_numbers <- rownames(nmrs_latlon)
# nmrs_latlon$temp_row <- result
# mean_temp_final2$row_numbers <- rownames(mean_temp_final2)
# 
# lat_lon_final <- merge(nmrs_latlon, mean_temp_final2, by.x="temp_row", by.y="row_numbers")
# lat_lon_final <- lat_lon_final[,-c(1,4:6)]
# colnames(lat_lon_final) <- c("lon_centre", "lat_centre", "summer_temperature")
# 
# ## map of mean air temp at NMRS 1km sqaures/points
# worldmap = map_data('world')
# ggplot() + 
#   geom_polygon(data = worldmap, 
#                aes(x = long, y = lat, group = group), 
#                fill = 'gray90', color = 'black') + 
#   coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#   geom_point(data = lat_lon_final, 
#              aes(x = as.numeric(lon_centre), 
#                  y = as.numeric(lat_centre), colour = summer_temperature), size=1) + 
#   scale_color_viridis_c(name="Mean summer temperature 1975-1991") + 
#   theme_void() +
#   theme(title = element_text(size = 12))
# 
# ## merge in with NMRS data
# nmrsdata_temp <- merge(nmrsdata, lat_lon_final, by=c("lon_centre", "lat_centre"))
# 
# ## save temp data
# saveRDS(nmrsdata_temp, file="../../Data/NMRS_data_summer_temperature.rds")
# 
# ##
# 
# ## now plot each species range of temperature across all occupied hectads
# ## min/max/median
# 
# ## need to subset by year (temperature values only refer to 1975-1991 so occupied hectads need to reflect this too)
# ## maybe should've been done earlier? Doesn't make any difference though
# 
# nmrsdata_temp <- nmrsdata_temp[nmrsdata_temp$Year>=1975 & nmrsdata_temp$Year<=1991, ]
# unique(nmrsdata_temp$Year) ## 17 years
# 
# detach(package:plyr)
# library(dplyr)
# library(ggplot2)
# summary_nmrs_temp <- nmrsdata_temp %>% 
#   group_by(Common_Name) %>% 
#   summarise(lower = min(summer_temperature), upper = max(summer_temperature), p = median(summer_temperature)) ## 796 species
# summary_nmrs_temp$lower <- as.numeric(summary_nmrs_temp$lower)
# summary_nmrs_temp$upper <- as.numeric(summary_nmrs_temp$upper)
# summary_nmrs_temp$p <- as.numeric(summary_nmrs_temp$p)
# 
# summary_nmrs_temp <- summary_nmrs_temp[order(summary_nmrs_temp$p, decreasing = TRUE),]  
# library(Hmisc)
# summary_nmrs_temp$groups<-as.numeric(cut2(summary_nmrs_temp$p, g=12))
# ## now make one plot for each group so they are readable 
# group <- unique(summary_nmrs_temp$groups)
# 
# for(i in group) {
#   print(i)
#   # Printing ggplot within for-loop
#   
#   temp_plot <- ggplot(data = summary_nmrs_temp[summary_nmrs_temp$groups==i,], mapping = aes(x = reorder(Common_Name, -p), y = p)) +
#     geom_pointrange(mapping = aes(ymin = lower, ymax = upper)) +
#     coord_flip() +
#     labs(x="Common name", y="Median summer temperature") +
#     ylim(9,18) +
#     theme_classic()
#   ggsave(temp_plot, file=paste0("../../Data/UKCP_tas_monthly/NMRS_summer_temp_range_", i,".png"), width = 20, height = 25, units = "cm")
#   Sys.sleep(2)
# }
# 
# library(plyr)
# 
# count(summary_nmrs_temp, "groups") ## 58-77 species in each group
# 
# ## all species - unreadable!
# ggplot(data = summary_nmrs_temp, mapping = aes(x = reorder(Common_Name, -p), y = p)) +
#   geom_pointrange(mapping = aes(ymin = lower, ymax = upper)) +
#   coord_flip() +
#   labs(x="Common name", y="Median summer temperature") +
#   theme_classic()
# #
# 
# 
# ###########################################################################################################################
# 
# ###############
# ## 1km scale ##
# ###############
# 
# ## Annual temperature
# library(ncdf4)
# library(raster)
# files <- list.files(path="../../Data/UKCP_tas_annual_1km/",
#                     pattern='*.nc',full.names=TRUE)
# 
# s <- raster::stack(files)
# #s[s %in% 9.969209968386869047442886268468442020e+36] <- NA ## change these values to NA
# 
# mean <- stackApply(s, indices =  rep(1,nlayers(s)), fun = "mean", na.rm = T) ## take mean across the 17 years/files
# plot(mean) # looks good
# 
# mean_newproj <- projectRaster(mean,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ## reproject to correct projection
# plot(mean_newproj, axes = T) ## mean temperature averaged across years 1975-1991 in 10x10km squares 
# 
# ##
# 
# mean_temp <- as.data.frame(rasterToPoints(mean_newproj)) ## get temp values as dataframe (246882 rows)
# colnames(mean_temp) <- c("lon", "lat", "temperature")
# 
# ## better plot of temperature data
# library(ggplot2)
# worldmap = map_data('world')
# temp_map <- ggplot() + 
#   geom_polygon(data = worldmap, 
#                aes(x = long, y = lat, group = group), 
#                fill = 'gray90', color = 'black') + 
#   coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#   geom_point(data = mean_temp, 
#              aes(x = as.numeric(lon), 
#                  y = as.numeric(lat), colour = temperature), size=1) + 
#   scale_color_viridis_c(name="Mean air temperature") + 
#   theme_void() +
#   theme(title = element_text(size = 12))
# temp_map
# 
# ## now match these with NMRS data
# nmrsdata <- readRDS("../../Data/NMRS/NMRS_grid_square_cleaned.rds")
# ## add 500m to easting and northing values to move record to centre of grid square
# ## this will ensure we are taking the temperature which represents middle of grid square
# nmrsdata$easting_centre <- nmrsdata$easting + 500
# nmrsdata$northing_centre <- nmrsdata$northing + 500
# ## get lat/lon values from easting and northing centre values
# wgs84 = "+init=epsg:4326"
# bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 
# +ellps=airy +datum=OSGB36 +units=m +no_defs'
# 
# ConvertCoordinates <- function(easting,northing) {
#   out = cbind(easting,northing)
#   mask = !is.na(easting)
#   sp <-  sp::spTransform(sp::SpatialPoints(list(easting[mask],northing[mask]),proj4string=sp::CRS(bng)),sp::CRS(wgs84))
#   out[mask,]=sp@coords
#   out
# }
# 
# x <- ConvertCoordinates(nmrsdata$easting_centre, nmrsdata$northing_centre)
# colnames(x) <- c("lon_centre","lat_centre")
# nmrsdata <- cbind(nmrsdata, x)
# 
# ## now match centre lat/lon NMRS values with nearest mean_temp lat/lon values
# nmrs_latlon <- nmrsdata[,c(11:12)]
# nmrs_latlon <- data.frame(unique(nmrs_latlon))
# 
# ## need to source Find Closest Rcpp function first!
# 
# find_closest <- function(lat1, long1, lat2, long2) {
#   
#   toRad <- pi / 180
#   lat1  <- lat1  * toRad
#   long1 <- long1 * toRad
#   lat2  <- lat2  * toRad
#   long2 <- long2 * toRad
#   
#   ord1  <- order(lat1)
#   rank1 <- match(seq_along(lat1), ord1)
#   ord2  <- order(lat2)
#   
#   ind <- find_closest_point(lat1[ord1], long1[ord1], lat2[ord2], long2[ord2])
#   
#   ord2[ind + 1][rank1]
# }
# 
# 
# result <- find_closest(nmrs_latlon[, 2], nmrs_latlon[, 1], 
#                        mean_temp[, 2], mean_temp[, 1])
# 
# result ### works super fast!!! assume numbers are the corresponding row numbers that each lat_lon value matches in temp_latlon?
# ## 62539 results, but 59967 unique values - some temperature sites have matched with multiple NMRS sites (mainly on the coast & Spurn Head)
# ## add row numbers to temp_latlon
# nmrs_latlon$row_numbers <- rownames(nmrs_latlon)
# nmrs_latlon$temp_row <- result
# mean_temp$row_numbers <- rownames(mean_temp)
# 
# lat_lon_final <- merge(nmrs_latlon, mean_temp, by.x="temp_row", by.y="row_numbers")
# lat_lon_final <- lat_lon_final[,-c(1,4:6)]
# colnames(lat_lon_final) <- c("lon_centre", "lat_centre", "temperature")
# 
# 
# ## map of mean air temp at NMRS 1km sqaures/points
# worldmap = map_data('world')
# ggplot() + 
#   geom_polygon(data = worldmap, 
#                aes(x = long, y = lat, group = group), 
#                fill = 'gray90', color = 'black') + 
#   coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#   geom_point(data = lat_lon_final, 
#              aes(x = as.numeric(lon_centre), 
#                  y = as.numeric(lat_centre), colour = temperature), size=1) + 
#   scale_color_viridis_c(name="Mean temperature 1975-1991") + 
#   theme_void() +
#   theme(title = element_text(size = 12))
# 
# ## merge in with NMRS data
# nmrsdata_temp <- merge(nmrsdata, lat_lon_final, by=c("lon_centre", "lat_centre"))
# 
# ## save temp data
# saveRDS(nmrsdata_temp, file="../../Data/NMRS_data_temperature_1km.rds")
# nmrsdata_temp <- readRDS("../../Data/NMRS_data_temperature_1km.rds")
# ##
# 
# ## now plot each species range of temperature across all occupied hectads
# ## min/max/median
# 
# ## need to subset by year (temperature values only refer to 1975-1991 so occupied hectads need to reflect this too)
# ## maybe should've been done earlier? Doesn't make any difference though
# 
# nmrsdata_temp <- nmrsdata_temp[nmrsdata_temp$Year>=1975 & nmrsdata_temp$Year<=1991, ]
# unique(nmrsdata_temp$Year) ## 17 years
# nmrsdata_temp$Year <- NULL
# nmrsdata_temp <- unique(nmrsdata_temp)
# 
# detach(package:plyr)
# library(dplyr)
# library(ggplot2)
# summary_nmrs_temp <- nmrsdata_temp %>% 
#   group_by(Common_Name) %>% 
#   summarise(lower = min(temperature), upper = max(temperature), p = median(temperature)) ## 793 species
# summary_nmrs_temp$lower <- as.numeric(summary_nmrs_temp$lower)
# summary_nmrs_temp$upper <- as.numeric(summary_nmrs_temp$upper)
# summary_nmrs_temp$p <- as.numeric(summary_nmrs_temp$p)
# 
# summary_nmrs_temp <- summary_nmrs_temp[order(summary_nmrs_temp$p, decreasing = TRUE),]  
# library(Hmisc)
# summary_nmrs_temp$groups<-as.numeric(cut2(summary_nmrs_temp$p, g=12))
# ## now make one plot for each group so they are readable 
# group <- unique(summary_nmrs_temp$groups)
# 
# sp_groups <- summary_nmrs_temp[,c(1,5)]
# 
# library(Hmisc)
# medians <- nmrsdata_temp %>% group_by(Common_Name) %>% summarise(median_temp = median(temperature))
# nmrsdata_temp <- merge(nmrsdata_temp, medians, by="Common_Name")
# ## now order by median
# nmrsdata_temp$median_temp <- as.numeric(as.character(nmrsdata_temp$median_temp))
# nmrsdata_temp <- nmrsdata_temp[order(nmrsdata_temp$median_temp),]
# ## now make one plot for each group so they are readable 
# nmrsdata_temp$group1 <- as.numeric(factor(nmrsdata_temp$Common_Name))
# nmrsdata_temp$group2<-as.numeric(cut2(nmrsdata_temp$group1, g=12))
# group <- unique(nmrsdata_temp$group2)
# 
# nmrsdata_temp <- merge(nmrsdata_temp, sp_groups, by="Common_Name", all=TRUE)
# group <- unique(nmrsdata_temp$groups)
# 
# sp1 <- nmrsdata[,c(3)]
# sp1 <- unique(sp1)
# sp2 <- 
#   indx <- match(nmrsdata$Common_Name, summary_nmrs_temp$Common_Name, nomatch = 0)
# 
# for(i in group) {
#   print(i)
#   # Printing ggplot within for-loop
#   
#   temp_plot <- ggplot(data = summary_nmrs_temp[summary_nmrs_temp$groups==i,], mapping = aes(x = reorder(Common_Name, -p), y = p)) +
#     geom_pointrange(mapping = aes(ymin = lower, ymax = upper)) +
#     coord_flip() +
#     labs(x="Common name", y="Median temperature") +
#     ylim(1,12) +
#     theme_classic()
#   # temp_plot <- ggplot(data = nmrsdata_temp[nmrsdata_temp$groups==i,], mapping = aes(x = reorder(Common_Name,-temperature, FUN=median), y = temperature)) +
#   #   geom_boxplot() +
#   #   coord_flip() +
#   #   labs(x="Common name", y="Median temperature") +
#   #   stat_summary(fun=mean, geom="point", shape=20, size=2, color="red", fill="red") +
#   #   ylim(3,11) +
#   #   theme_classic()
#   
#   ggsave(temp_plot, file=paste0("../../Data/UKCP_tas_annual_1km/NMRS_1km_temp_range_", i,".png"), width = 20, height = 25, units = "cm")
#   Sys.sleep(2)
# }
# 
# library(plyr)
# 
# count(summary_nmrs_temp, "groups") ## 65-67 species in each group
# 
# ## all species - unreadable!
# ggplot(data = summary_nmrs_temp, mapping = aes(x = reorder(Common_Name, -p), y = p)) +
#   geom_pointrange(mapping = aes(ymin = lower, ymax = upper)) +
#   coord_flip() +
#   labs(x="Common name", y="Median temperature") +
#   theme_classic()
#



