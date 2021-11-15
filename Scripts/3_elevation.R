##########################
#### user: Lisbeth Hordley
#### date: July 2021
#### info: Add in elevation data for NMRS sites

rm(list = ls())

## load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(sf)
library(ggplot2)
library(plyr)
library(ggpubr)


## elevation data
library(raster)
tmp = raster("Data/elevation1x1_new.tif")
crs(tmp) ## projection
res(tmp) ## resolution (1000 x 1000m) (needs multiplying by a factor of 10 to be 10,000 x 10,000m)
tmp_agg <- aggregate(tmp, fact=10, fun=mean)
res(tmp_agg) ## 10000 x 10000m = 10x10km
plot(tmp_agg)
hist(values(tmp_agg))

## re-project using NMRS cleaned hectad data
nmrsdata <- readRDS("Data/NMRS/NMRS_hectad_cleaned.rds") 
nmrs_lat_lon <- unique(nmrsdata[,c(1,7:8)]) # 2719 (keep gridref to merge elevation with nmrsdata later)

tmp_newproj <- projectRaster(tmp_agg,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(tmp_newproj, axes = T)

tmp_elev <- raster::extract(tmp_newproj, nmrs_lat_lon[,2:3], df=T)
tmp_elev <- cbind(tmp_elev, nmrs_lat_lon)
## multiply elevation by 10 - see metadata https://www.eea.europa.eu/data-and-maps/data/digital-elevation-model-of-europe?fbclid=IwAR1yqtoCIjiV1QIe30h8DCaYhWH-7GO0X68CKLAkD5FY9e5WMs6OPU9YSas 
tmp_elev$elevation1x1_new <- tmp_elev$elevation1x1_new*10 ## elevation in meters 
tmp_elev <- tmp_elev[,-1]
colnames(tmp_elev)[1] <- "elevation10x10km"

## map of elevation at NMRS 1km sqaures/points
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = tmp_elev, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour = elevation10x10km), size=1) + 
  scale_color_viridis_c(name="Elevation") + 
  theme_void() +
  theme(title = element_text(size = 12))



##############################################################
### Get standard deviation at hectad level using 1km elevation
elev_sd <- aggregate(tmp, fact=10, fun=sd)
res(elev_sd) ## 10000 x 10000m = 10x10km
plot(elev_sd)

elev_sd_newproj <- projectRaster(elev_sd,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(elev_sd_newproj, axes = T)

elev_sd_nmrs <- raster::extract(elev_sd_newproj, nmrs_lat_lon[,2:3], df=T)
elev_sd_nmrs <- cbind(elev_sd_nmrs, nmrs_lat_lon)
## multiply elevation by 10 - see metadata https://www.eea.europa.eu/data-and-maps/data/digital-elevation-model-of-europe?fbclid=IwAR1yqtoCIjiV1QIe30h8DCaYhWH-7GO0X68CKLAkD5FY9e5WMs6OPU9YSas 
elev_sd_nmrs$elevation1x1_new <- elev_sd_nmrs$elevation1x1_new*10 ## elevation in meters 
elev_sd_nmrs <- elev_sd_nmrs[,-1]
colnames(elev_sd_nmrs)[1] <- "elevation10x10km_SD"

## plot standard deviation of elevation
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = elev_sd_nmrs, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour = elevation10x10km_SD), size=1) + 
  scale_color_viridis_c(name="SD of elevation") + 
  theme_void() +
  theme(title = element_text(size = 12))

###

## remove lat/lon - not needed and doesn't merge in with nmrs data properly
tmp_elev <- tmp_elev[,-c(3:4)]
elev_sd_nmrs <- elev_sd_nmrs[,-c(3:4)]
## put elevation & standard deviation together
elev_mean_sd <- merge(tmp_elev, elev_sd_nmrs, by="Hectad")
## if throwing memory issues: memory.limit(size = 15000)
nmrsdata <- merge(nmrsdata, elev_mean_sd, by="Hectad", all=T)
## 20,501,968 rows in NMRS data (same as when first read in)
head(nmrsdata)
length(unique(nmrsdata$Hectad)) ## 2719
## save data
saveRDS(nmrsdata, file="Data/NMRS/NMRS_hectad_elevation.rds")
## use this for trailing edge analysis for elevation and latitude at a 10km scale

hist(elev_sd_nmrs$elevation10x10km_SD)
## exclude the tail of standard deviation - this will remove hectads with large variance in elevation
## therefore those hectads for which the elevation value is less likely to represent the actual elevation
## where the species is found

## which species are found in hectads with high variance? 
