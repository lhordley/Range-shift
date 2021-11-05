##########################
#### user: Lisbeth Hordley
#### date: July 2021
#### info: Adding in elevation data for NMRS sites

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
tmp = raster("../../Range shift/Data/elevation1x1_new.tif")
crs(tmp) ## projection
res(tmp) ## resolution (1000 x 1000m) (needs multiplying by a factor of 10 to be 10,000 x 10,000m)
tmp_agg <- aggregate(tmp, fact=10)
res(tmp_agg) ## 10000 x 10000m = 10x10km
plot(tmp_agg)
hist(values(tmp_agg))

## re-project using NMRS cleaned hectad data
nmrsdata <- readRDS("../../Range shift/Data/NMRS/NMRS_hectad_cleaned.rds") 
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

## remove lat/lon - not needed and doesn't merge in with nmrs data properly
tmp_elev <- tmp_elev[,-c(3:4)]
## if throwing memory issues: memory.limit(size = 15000)
nmrsdata <- merge(nmrsdata, tmp_elev, by="Hectad", all=T)
## 20,501,968 rows in NMRS data (same as when first read in)
head(nmrsdata)
length(unique(nmrsdata$Hectad))
## save data
saveRDS(nmrsdata, file="../../Range shift/Data/NMRS/NMRS_hectad_elevation.rds")
## use this for trailing edge analysis for elevation and latitude at a 10km scale

