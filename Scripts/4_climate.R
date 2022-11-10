##########################
#### user: Lisbeth Hordley
#### date: October 2021
#### info: Add climate data for NMRS sites

rm(list = ls())

## load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(sf)
library(ggplot2)
library(plyr)
library(here)
library(raster)
library(gdata)

##############################################################################
########### MEAN ANNUAL TEMPERATURE EARLY TIME PERIOD: 1975 - 1991 ###########
##############################################################################

## Extract annual temperature data & merge with NMRS hectads

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

#############################################################################
########### MEAN ANNUAL TEMPERATURE LATE TIME PERIOD: 2012 - 2016 ###########
#############################################################################

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


##############################################################################
########### MEAN SUMMER TEMPERATURE EARLY TIME PERIOD: 1975 - 1991 ###########
##############################################################################

## read in temperature data (5km scale from metoffice CP18, monthly temperature for years 1975:1991)
files <- list.files(path="Data/UKCP_tas_monthly/Early_TP_75_91/",
                    pattern='*.nc',full.names=TRUE)

s <- raster::stack(files)
#s[s %in% 9.969209968386869047442886268468442020e+36] <- NA ## change these values to NA

s2 = raster(extent(s), resolution = 10000, crs = crs(s)) ## change resolution to 10x10km
s3 = raster::resample(s, s2)

mean_newproj <- projectRaster(s3,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ## reproject to correct projection
plot(mean_newproj, axes = T) ## monthly temp for each of the 17 years (haven't taken mean yet)

## June, July & August = Summer months
mean_temp <- as.data.frame(rasterToPoints(mean_newproj)) ## get temp values as dataframe (2929 rows, 206 columns (2 of them are lat/lon))
## need to rep years (1975:1991) and rep months (x12) to make this into a long dataframe
library(reshape2)
mean_temp_final <- melt(data=mean_temp, id.vars = c("x", "y"), variable.name = "dataset", value.name = "mean_temp")
months = rep(c("Jan","Feb","March","April","May","Jun","Jul","Aug","Sept","Oct",
               "Nov", "Dec"), times=17, each=2929)
years <- rep(1975:1991, each=35148)
mean_temp_final$month <- months
mean_temp_final$year <- years
## remove dataset
mean_temp_final$dataset <- NULL
colnames(mean_temp_final)[1:2] <- c("lon","lat")

## now only keep summer months
mean_temp_final <- subset(mean_temp_final, (month %in% c("Jun","Jul","Aug")))
## now take mean across summer months and across years (left with one temperature value for each hectad)
library(dplyr)
mean_temp_final2 <- mean_temp_final %>% group_by(lon, lat) %>% dplyr::summarize(summer_mean_temp=mean(mean_temp))
mean_temp_final2 <- as.data.frame(mean_temp_final2)

## save file
write.csv(mean_temp_final2, file="Data/UKCP_tas_monthly/Early_TP_75_91/mean_summer_temp_75_91.csv", row.names=FALSE)

## better plot of temperature data
library(ggplot2)
worldmap = map_data('world')
summer_temp_map <- ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  geom_point(data = mean_temp_final2,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat), colour = summer_mean_temp), size=1) +
  scale_color_viridis_c(name="Mean summer temperature") +
  theme_void() +
  theme(title = element_text(size = 12))
summer_temp_map

#############################################################################
########### MEAN SUMMER TEMPERATURE LATE TIME PERIOD: 2012 - 2016 ###########
#############################################################################

## read in temperature data (5km scale from metoffice CP18, monthly temperature for years 1975:1991)
files <- list.files(path="Data/UKCP_tas_monthly/Late_TP_12_16/",
                    pattern='*.nc',full.names=TRUE)

s <- raster::stack(files)
#s[s %in% 9.969209968386869047442886268468442020e+36] <- NA ## change these values to NA

s2 = raster(extent(s), resolution = 10000, crs = crs(s)) ## change resolution to 10x10km
s3 = raster::resample(s, s2)

mean_newproj <- projectRaster(s3,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ## reproject to correct projection
plot(mean_newproj, axes = T) ## monthly temp for each of the 5 years (haven't taken mean yet)

## June, July & August = Summer months
mean_temp <- as.data.frame(rasterToPoints(mean_newproj)) ## get temp values as dataframe (2929 rows, 62 columns (2 of them are lat/lon))
## need to rep years (2012:2016) and rep months (x12) to make this into a long dataframe
library(reshape2)
mean_temp_final <- melt(data=mean_temp, id.vars = c("x", "y"), variable.name = "dataset", value.name = "mean_temp")
months = rep(c("Jan","Feb","March","April","May","Jun","Jul","Aug","Sept","Oct",
               "Nov", "Dec"), times=5, each=2929)
years <- rep(2012:2016, each=35148)
mean_temp_final$month <- months
mean_temp_final$year <- years
## remove dataset
mean_temp_final$dataset <- NULL
colnames(mean_temp_final)[1:2] <- c("lon","lat")

## now only keep summer months
mean_temp_final <- subset(mean_temp_final, (month %in% c("Jun","Jul","Aug")))
## now take mean across summer months and across years (left with one temperature value for each hectad)
library(dplyr)
mean_temp_final2 <- mean_temp_final %>% group_by(lon, lat) %>% dplyr::summarize(summer_mean_temp=mean(mean_temp))
mean_temp_final2 <- as.data.frame(mean_temp_final2)

## save file
write.csv(mean_temp_final2, file="Data/UKCP_tas_monthly/Late_TP_12_16/mean_summer_temp_12_16.csv", row.names=FALSE)

## better plot of temperature data
library(ggplot2)
worldmap = map_data('world')
summer_temp_map <- ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  geom_point(data = mean_temp_final2,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat), colour = summer_mean_temp), size=1) +
  scale_color_viridis_c(name="Mean summer temperature") +
  theme_void() +
  theme(title = element_text(size = 12))
summer_temp_map



######################################################################################
########### MEAN ANNUAL TOTAL PRECIPITATION EARLY TIME PERIOD: 1975 - 1991 ###########
######################################################################################

## Extract annual precipitation data & merge with NMRS hectads

## read in precipitation data (5km scale from metoffice CP18, total precipitation for years 1975:1991)
library(ncdf4)
library(raster)
files <- list.files(path="Data/UKCP_precip_annual/Early_TP_75_91/",
                    pattern='*.nc',full.names=TRUE)

s <- raster::stack(files)

s2 = raster(extent(s), resolution = 10000, crs = crs(s)) ## change resolution to 10x10km
s3 = raster::resample(s, s2)

mean <- stackApply(s3, indices =  rep(1,nlayers(s3)), fun = "mean", na.rm = T) ## take mean across the 17 years/files
plot(mean) # looks good

mean_newproj <- projectRaster(mean,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ## reproject to correct projection
plot(mean_newproj, axes = T) ## mean precipitation averaged across years 1975-1991 in 10x10km squares 

##

mean_precip <- as.data.frame(rasterToPoints(mean_newproj)) ## get temp values as dataframe (2929 rows)
colnames(mean_precip) <- c("lon", "lat", "total_precip")

## save file
write.csv(mean_precip, file="Data/UKCP_precip_annual/Early_TP_75_91/mean_total_precip_75_91.csv", row.names=FALSE)

## better plot of precipitation data
library(ggplot2)
worldmap = map_data('world')
precip_map <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = mean_precip, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour = total_precip), size=1) + 
  scale_color_viridis_c(name="Mean total precipitation") + 
  theme_void() +
  theme(title = element_text(size = 12))
precip_map


#####################################################################################
########### MEAN ANNUAL TOTAL PRECIPITATION LATE TIME PERIOD: 2012 - 2016 ###########
#####################################################################################

## read in precipitation data (5km scale from metoffice CP18, total precipitation for years 1975:1991)
library(ncdf4)
library(raster)
files <- list.files(path="Data/UKCP_precip_annual/Late_TP_12_16/",
                    pattern='*.nc',full.names=TRUE)

s <- raster::stack(files)

s2 = raster(extent(s), resolution = 10000, crs = crs(s)) ## change resolution to 10x10km
s3 = raster::resample(s, s2)

mean <- stackApply(s3, indices =  rep(1,nlayers(s3)), fun = "mean", na.rm = T) ## take mean across the 17 years/files
plot(mean) # looks good

mean_newproj <- projectRaster(mean,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") ## reproject to correct projection
plot(mean_newproj, axes = T) ## mean precipitation averaged across years 1975-1991 in 10x10km squares 

##

mean_precip <- as.data.frame(rasterToPoints(mean_newproj)) ## get temp values as dataframe (2929 rows)
colnames(mean_precip) <- c("lon", "lat", "total_precip")

## save file
write.csv(mean_precip, file="Data/UKCP_precip_annual/Late_TP_12_16/mean_total_precip_12_16.csv", row.names=FALSE)

## better plot of precipitation data
library(ggplot2)
worldmap = map_data('world')
precip_map <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = mean_precip, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour = total_precip), size=1) + 
  scale_color_viridis_c(name="Mean total precipitation") + 
  theme_void() +
  theme(title = element_text(size = 12))
precip_map


##################################################################################################
##################################################################################################
##################################################################################################
rm(list = ls())
options(scipen=999)


## Put climate datasets together and match with closest NMRS sites
mean_temp_tp1 <- read.csv("Data/UKCP_tas_annual/Early_TP_75_91/mean_temp_75_91.csv", header=TRUE)
mean_temp_tp2 <- read.csv("Data/UKCP_tas_annual/Late_TP_12_16/mean_temp_12_16.csv", header=TRUE)
summer_temp_tp1 <- read.csv("Data/UKCP_tas_monthly/Early_TP_75_91/mean_summer_temp_75_91.csv", header=TRUE, strip.white = TRUE)
summer_temp_tp2 <- read.csv("Data/UKCP_tas_monthly/Late_TP_12_16/mean_summer_temp_12_16.csv", header=TRUE, strip.white = TRUE)
total_precip_tp1 <- read.csv("Data/UKCP_precip_annual/Early_TP_75_91/mean_total_precip_75_91.csv", header=TRUE)
total_precip_tp2 <- read.csv("Data/UKCP_precip_annual/Late_TP_12_16/mean_total_precip_12_16.csv", header=TRUE)

mean_temp_tp1$Time_period <- "1975-1991"
mean_temp_tp2$Time_period <- "2012-2016"
mean_temp <- rbind(mean_temp_tp1, mean_temp_tp2)
summer_temp_tp1$Time_period <- "1975-1991"
summer_temp_tp2$Time_period <- "2012-2016"
summer_temp <- rbind(summer_temp_tp1, summer_temp_tp2)
total_precip_tp1$Time_period <- "1975-1991"
total_precip_tp2$Time_period <- "2012-2016"
total_precip <- rbind(total_precip_tp1, total_precip_tp2)

mean_temp <- mean_temp[order(mean_temp$lat),]
summer_temp <- summer_temp[order(summer_temp$lat),]
total_precip <- total_precip[order(total_precip$lat),]
## need to round lon values - for some reason summer temp ones are 1 digit different 
mean_temp$lon <- round(mean_temp$lon ,8)
summer_temp$lon <- round(summer_temp$lon ,8)
total_precip$lon <- round(total_precip$lon ,8)

## put data together
climate_final <- merge(mean_temp, total_precip, by.x=c("lon", "lat", "Time_period"), 
                       by.y=c("lon", "lat", "Time_period"), all=TRUE)

climate_final <- merge(climate_final, summer_temp, by.x=c("lon", "lat", "Time_period"), 
                       by.y=c("lon", "lat", "Time_period"), all=TRUE)

## save file
write.csv(climate_final, file="Data/Climate_final.csv", row.names=FALSE)


## now match these with NMRS data
nmrsdata <- readRDS("Data/NMRS/NMRS_hectad_elevation.rds")
climate_final <- read.csv("Data/Climate_final.csv", header=TRUE)
## now match centre lat/lon NMRS values with nearest mean_temp lat/lon values


nmrs_latlon <- nmrsdata[,c("lat_centre", "lon_centre")] ## take lat and lon values from NMRS data
nmrs_latlon <- data.frame(unique(nmrs_latlon)) ## unique the dataframe to give one row per site (2695 rows)
climate_latlon <- climate_final[,c("lat", "lon")] ## same for climate data
climate_latlon <- data.frame(unique(climate_latlon)) ## 2695 rows
row.names(climate_latlon) <- NULL 
row.names(climate_latlon) <-1:nrow(climate_latlon) ## this makes sure the row numbers are sequential after taking unique of the df

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


result <- find_closest(nmrs_latlon[,"lat_centre"], nmrs_latlon[, "lon_centre"], 
                       climate_latlon[, "lat"], climate_latlon[, "lon"]) ## give the lat/lon values

# result ### works super fast!!! assume numbers are the corresponding row numbers that each lat_lon value matches in temp_latlon?
## 62539 results, but 59967 unique values - some temperature sites have matched with multiple NMRS sites (mainly on the coast & Spurn Head)
## add row numbers to temp_latlon
nmrs_latlon$row_numbers <- rownames(nmrs_latlon)
nmrs_latlon$climate_row <- result
climate_latlon$row_numbers <- rownames(climate_latlon)
climate_latlon$row_numbers <- as.integer(climate_latlon$row_numbers)

lat_lon_final <- merge(nmrs_latlon, climate_latlon, by.x="climate_row", by.y="row_numbers", all.x=TRUE)
## add in climate data and time period info again
lat_lon_final <- merge(lat_lon_final, climate_final, by=c("lat", "lon"))

lat_lon_final <- lat_lon_final[,-c(1:3,6)]

hectad <- unique(nmrsdata[,c("lat_centre", "lon_centre", "Hectad")])
lat_lon_final <- merge(lat_lon_final, hectad, by=c("lat_centre", "lon_centre"))
lat_lon_final <- lat_lon_final[,-c(1:2)]

## merge in with NMRS data
nmrsdata_climate <- merge(nmrsdata, lat_lon_final, by=c("Hectad", "Time_period"), all.x=TRUE)

saveRDS(nmrsdata_climate, file="Data/NMRS/NMRS_hectad_elevation_climate.rds") ## 667521 rows




########################################################################################################
## PLOT SPECIES TEMPERATURE RANGES (MEDIAN & 75TH PERCENTILE) BETWEEN 1975-1991
## These are given to Richard to determine which cool-adapted species to select
rm(list = ls())

nmrsdata_climate <- readRDS("Data/NMRS/NMRS_hectad_elevation_climate.rds") ## NMRS data with temperature 
migrant_hectads <- read.csv("Data/NMRS/NMRS_migrant_hectads_exclude.csv", header=TRUE) ## species x hectad combinations to remove due to immigrant populations
## provided by Richard in Teams 20/10/21

## Remove TP1 data - only looking at temperature range between 1975 - 1991
nmrsdata_climate <- nmrsdata_climate[nmrsdata_climate$Time_period=="1975-1991",] # 254689
migrant_hectads <- migrant_hectads[migrant_hectads$Time_period=="1975-1991",]
## merge migrant exclusion data with nmrs_temp data
## remove hectad/species combinations where exclude == 1
migrant_hectads[is.na(migrant_hectads)] <- 0
nmrsdata_climate <- merge(nmrsdata_climate, migrant_hectads, by=c("Common_name", "Hectad", "Time_period"), all=TRUE)
nmrsdata_climate[is.na(nmrsdata_climate)] <- 0
nmrsdata_climate <- nmrsdata_climate[nmrsdata_climate$Exclude !=1, ] ## 254544 rows

detach(package:plyr)
library(dplyr)
library(ggplot2)
summary_nmrs_temp <- nmrsdata_climate %>% 
  group_by(Common_name) %>% 
  summarise(lower = min(temperature), upper = max(temperature), p = median(temperature)) ## 788 species
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
nmrsdata_climate <- merge(nmrsdata_climate, sp_groups, by="Common_name", all=TRUE)
group <- unique(nmrsdata_climate$groups) ## 12 groups

for(i in group) {
  print(i)
  # Printing ggplot within for-loop
  
  # temp_plot <- ggplot(data = summary_nmrs_temp[summary_nmrs_temp$groups==i,], mapping = aes(x = reorder(Common_Name, -p), y = p)) +
  #   geom_pointrange(mapping = aes(ymin = lower, ymax = upper)) +
  #   coord_flip() +
  #   labs(x="Common name", y="Median temperature") +
  #   ylim(3,12) +
  #   theme_classic()
  temp_plot <- ggplot(data = nmrsdata_climate[nmrsdata_climate$groups==i,], mapping = aes(x = reorder(Common_name,-temperature, FUN=median), y = temperature)) +
    geom_boxplot() +
    #geom_hline(yintercept=9, linetype="dotted") +
    coord_flip() +
    labs(x="Common name", y="Median temperature") +
    scale_y_continuous(breaks = seq(3, 11, by = 1)) +
    theme_light()
  
  ggsave(temp_plot, file=paste0("Data/UKCP_tas_annual/Early_TP_75_91/NMRS_new_temp_range_box_", i,".png"), width = 20, height = 25, units = "cm")
  Sys.sleep(2)
}

## calculate median and 75th percentile for each species for Richard
temp_ranges <- nmrsdata_climate %>%  group_by(Common_name) %>% 
  summarise(median_temperature = median(temperature), percentile_75th = quantile(temperature, probs = 0.75))
## order by median ascending
temp_ranges <- arrange(temp_ranges, median_temperature)
## save file
write.csv(temp_ranges, file="Data/NMRS_new_temp_median_percentile.csv", row.names=FALSE)
## Boxplots and percentile list given to Richard to select cool adapted species

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

