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
## add 5km to easting and northing values to move record to centre of hectad
## this will ensure we are taking the mean elevation which represents the mean of surrounding 1km squares WITHIN the hectad
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

nmrs_lat_lon <- unique(nmrsdata[,c(1,12,13)]) # 2695 (keep gridref to merge elevation with nmrsdata later)

tmp_newproj <- projectRaster(tmp_agg,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(tmp_newproj, axes = T)

tmp_elev <- raster::extract(tmp_newproj, nmrs_lat_lon[,2:3], df=T)
tmp_elev <- cbind(tmp_elev, nmrs_lat_lon)
## multiply elevation by 10 - see metadata https://www.eea.europa.eu/data-and-maps/data/digital-elevation-model-of-europe?fbclid=IwAR1yqtoCIjiV1QIe30h8DCaYhWH-7GO0X68CKLAkD5FY9e5WMs6OPU9YSas 
tmp_elev$elevation1x1_new <- tmp_elev$elevation1x1_new*10 ## elevation in meters 
tmp_elev <- tmp_elev[,-1]
colnames(tmp_elev)[1] <- "elevation10x10km"

## map of elevation at NMRS 10km sqaures/points
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = tmp_elev, 
             aes(x = as.numeric(lon_centre), 
                 y = as.numeric(lat_centre), colour = elevation10x10km), size=1) + 
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
elev_variation <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = elev_sd_nmrs, 
             aes(x = as.numeric(lon_centre), 
                 y = as.numeric(lat_centre), colour = elevation10x10km_SD), size=1) + 
  scale_color_viridis_c(name="SD of elevation") + 
  theme_void() +
  theme(title = element_text(size = 12))
elev_variation
ggsave(elev_variation, file="Maps/Elevation10x10km_SD.png")
###

## remove lat/lon - not needed and doesn't merge in with nmrs data properly
tmp_elev <- tmp_elev[,-c(3:4)]
elev_sd_nmrs <- elev_sd_nmrs[,-c(3:4)]
## put elevation & standard deviation together
elev_mean_sd <- merge(tmp_elev, elev_sd_nmrs, by="Hectad")
## if throwing memory issues: memory.limit(size = 15000)
nmrsdata <- merge(nmrsdata, elev_mean_sd, by="Hectad", all=T)
## 667,521 rows in NMRS data (same as when first read in)
head(nmrsdata)
length(unique(nmrsdata$Hectad)) ## 2695
## save data
saveRDS(nmrsdata, file="Data/NMRS/NMRS_hectad_elevation.rds")
## use this for trailing edge analysis for elevation and latitude at a 10km scale

hist(elev_sd_nmrs$elevation10x10km_SD)
## exclude the tail of standard deviation - this will remove hectads with large variance in elevation
## therefore those hectads for which the elevation value is less likely to represent the actual elevation
## where the species is found

## which species are found in hectads with high variance? 
