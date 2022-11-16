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
total_precip_tp1 <- read.csv("Data/UKCP_precip_annual/Early_TP_75_91/mean_total_precip_75_91.csv", header=TRUE)
total_precip_tp2 <- read.csv("Data/UKCP_precip_annual/Late_TP_12_16/mean_total_precip_12_16.csv", header=TRUE)

mean_temp_tp1$Time_period <- "1975-1991"
mean_temp_tp2$Time_period <- "2012-2016"
mean_temp <- rbind(mean_temp_tp1, mean_temp_tp2)
total_precip_tp1$Time_period <- "1975-1991"
total_precip_tp2$Time_period <- "2012-2016"
total_precip <- rbind(total_precip_tp1, total_precip_tp2)

mean_temp <- mean_temp[order(mean_temp$lat),]
total_precip <- total_precip[order(total_precip$lat),]
## need to round lon values - for some reason summer temp ones are 1 digit different 
mean_temp$lon <- round(mean_temp$lon ,8)
total_precip$lon <- round(total_precip$lon ,8)

## put data together
climate_final <- merge(mean_temp, total_precip, by.x=c("lon", "lat", "Time_period"), 
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




