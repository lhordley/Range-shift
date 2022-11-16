##########################
#### user: Lisbeth Hordley
#### date: October 2021
#### info: Multidimensional range shift analysis using NMRS data

rm(list = ls())
options(scipen=999)

## load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(sf)
library(ggplot2)
library(plyr)
library(ggpubr)
library(Rmisc)

## NMRS data for cool-adapted moths with elevation and temperature data at a 10km scale
nmrsdata <- readRDS("Data/NMRS/NMRS_cool_moths_final.rds") ## NMRS data for all hectads and all years with elevation
## hectad recording levels
hec_records <- read.csv("Data/NMRS/Hectad_recording_levels_1975_1991_2012_2016.csv", header=TRUE)

nmrsdata <- nmrsdata[!nmrsdata$Common_name=="Lunar Thorn",]

###################################
## well & heavily recorded hectads

## filter to recorded hectads only (i.e. those recorded once in both time periods)
well_heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded" | RECORDING.LEVEL=="Well recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
well_heavy_hecs <- well_heavy_hecs[which(well_heavy_hecs$n_row > 1), ]
well_heavy_hecs$n_row <- NULL ## 1084 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_well <- nmrsdata[which(nmrsdata$Hectad %in% well_heavy_hecs$HECTAD), ]
length(unique(nmrsdata_well$Hectad)) ## 962 hectads
length(unique(nmrsdata_well$Common_name)) ## 76 species

#### 1. Multidimensional RANGE shift 
# Distance and direction of centroid shift between TP1 and TP2
# take mean easting and northing and lat and lon across each species occupied hectads for each time period = range centroid
# following Gillings et al 2015
range_centroid_well <- nmrsdata_well %>% group_by(Common_name, Time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
range_centroid_well2 <- range_centroid_well %>%
  gather(key, value, -Common_name, -Time_period) %>%
  unite(col, key, Time_period) %>%
  spread(col, value)
## remove NAs - some species are only found in one TP - Grey and Scotch Burnet removed
range_centroid_well2 <- na.omit(range_centroid_well2) ## 74 species
# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
range_centroid_well2$distance <- sqrt((range_centroid_well2$`mean_northing_1975-1991`-range_centroid_well2$`mean_northing_2012-2016`)^2 + 
                                       (range_centroid_well2$`mean_easting_1975-1991`-range_centroid_well2$`mean_easting_2012-2016`)^2) ## pythagoras equation
range_centroid_well2$distance <- range_centroid_well2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
range_centroid_well2$bearing <- bearing(range_centroid_well2[,c(6,4)], range_centroid_well2[,c(7,5)])
range_centroid_well2$direction <- (range_centroid_well2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360

# mean_ci <- data.frame(Common_Name=c("mean", "upper", "lower"), distance=c(200,200,200), direction=c(320.1, 288.77, 339.76))
# 
# range_centroid_well2 <- rbind(range_centroid_well2, mean_ci)

## plot result
ybreaks <- c(0,50,100,150,200,250,300)
well_centroid <- ggplot(data=range_centroid_well2,
                        aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1), lwd=0.6) +
  #geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F) +
  theme_bw() +
  theme(text = element_text(size=14), panel.border = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
well_centroid
ggsave(well_centroid, file="Graphs/multidirectional_shifts_well_new.png")

## circular plot for poster
ybreaks <- c(0,50,100,150,200,250,300)
well_centroid2 <- ggplot(data=range_centroid_well2,
                        aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1), lwd=0.6) +
  geom_point(aes(x=322, y=320), colour="red", lwd=1) +
  #geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F) +
  theme_bw()
  #theme(text = element_text(size=14), panel.border = element_blank(), axis.text.y = element_blank(), 
        #axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
well_centroid2
ggsave(well_centroid2, file="Graphs/poster_result2.png")


worldmap = map_data('world')
worldmap <- worldmap[!worldmap$region=="Ireland",]
worldmap <- worldmap[!worldmap$subregion=="Northern Ireland",]

rh1 <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(aes(x=-2.2, y=54.2), size=4, colour="#D55E00") +
  geom_point(aes(x=-3.6, y=55.5), size=4, colour="#0072B2") +
  #geom_segment(aes(x = -2.4,y = 54.4,xend = -3.4,yend = 55.3),arrow=arrow(), lwd=1.5) +
  theme(title = element_text(size = 12))
rh1
ggsave(rh1, file="Graphs/poster_method2.png", width=15, height=14, units="cm")

ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_well[nmrsdata_well$Common_name=="Welsh Clearwing" & nmrsdata_well$Time_period=="1975-1991",], 
             aes(x = lon, y=lat, colour=Time_period), size=1) +
  
  theme(title = element_text(size = 12))


ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =range_centroid_well[range_centroid_well$Common_Name=="(Ruddy Highflyer)",], 
             aes(x = mean_lon, y=mean_lat, colour=time_period), size=1) +
  
  theme(title = element_text(size = 12))

# use circ.mean from the circstats package to find a mean direction & confidence limits across all species
library(NISTunits)
library(CircStats)
range_centroid_well2$radians <- NISTdegTOradian(range_centroid_well2$direction)
x <- deg(circ.mean(range_centroid_well2$radians)) # -36, add full circle
(x+360)%%360 ## 322
ci <- vm.bootstrap.ci(range_centroid_well2$radians) 
deg(c(-1.21, -0.27)) 
((-69.32789)+360)%%360 ## 291
((-15.46986)+360)%%360 ## 345

# median & quartiles for distances
median(range_centroid_well2$distance) ## 36.9km 
quantile(range_centroid_well2$distance) ## 25th = 17.75km, 75th = 55.95km

## Rayleigh's test
## if significant, the distribution of directions differ significantly from an even distribution
## if R = 0 (the data are completely spread around the circle)
## if R = 1 (the data are completely concentrated on one point)
r.test(range_centroid_well2$radians)
## p < 0.001
## r = 0.33
## low R = variable directions, but significant means it has significant bias to a certain direction (the mean)

#### 2. Multidimensional CLIMATE shift
## What direction should species move in to remain in the same climatic envelope? 
# Then match this with the direction and distance species have actually shifted
# Do they correlate? 

##############
# 2a: Mean annual temperature
##############

#### Direction and distance of shift based on mean annual temperature
## 2.5th and 97.5th percentiles for each species in TP1
detach(package:Rmisc)
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_well %>% group_by(Common_name) %>% filter(Time_period=="1975-1991") %>%
  dplyr::summarise(temp_2.5_perc = quantile(temperature, probs = 0.025),
            temp_97.5_perc = quantile(temperature, probs = 0.975))
nmrsdata_well <- merge(nmrsdata_well, tp1_species_percentiles, by="Common_name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_well %>% group_by(Common_name, Time_period) %>%
  filter(temperature >= temp_2.5_perc & temperature <= temp_97.5_perc)

## calculate centroids of each
temp_centroid_well <- filtered %>% group_by(Common_name, Time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing)) # 74 species

## change orientation of data
temp_centroid_well2 <- temp_centroid_well %>%
  gather(key, value, -Common_name, -Time_period) %>%
  unite(col, key, Time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2

## remove NAs
temp_centroid_well2 <- na.omit(temp_centroid_well2) ## 66 species (7 species removed)

sp_all <- data.frame(unique(nmrsdata_well$Common_name))
sp_sub <- data.frame(unique(temp_centroid_well2$Common_name))
sp_list <- sp_all %>% 
  filter(!sp_all$unique.nmrsdata_well.Common_name. %in% sp_sub$unique.temp_centroid_well2.Common_name.)
# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
temp_centroid_well2$distance <- sqrt((temp_centroid_well2$`mean_northing_1975-1991`-temp_centroid_well2$`mean_northing_2012-2016`)^2 + 
                                      (temp_centroid_well2$`mean_easting_1975-1991`-temp_centroid_well2$`mean_easting_2012-2016`)^2) ## pythagoras equation
temp_centroid_well2$distance <- temp_centroid_well2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
temp_centroid_well2$bearing <- bearing(temp_centroid_well2[,c(6,4)], temp_centroid_well2[,c(7,5)])
temp_centroid_well2$direction <- (temp_centroid_well2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
temp_centroid_well2$radians <- NISTdegTOradian(temp_centroid_well2$direction)
x <- deg(circ.mean(temp_centroid_well2$radians)) # -11, add full circle
(x+360)%%360 ## 341
ci <- vm.bootstrap.ci(temp_centroid_well2$radians) 
deg(c(5.79, 6.11)) ## 332 350

# median & quartiles for distances
median(temp_centroid_well2$distance) ## 77.19km
quantile(temp_centroid_well2$distance) ## 25th = 29.43km, 75th = 103.01km

## plot result
ybreaks <- c(0,50,100,150,200,250,300)
well_temp_centroid <- ggplot(data=temp_centroid_well2,
                             aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1), lwd=0.6) +
  #geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F) +
  theme_bw() +
  theme(text = element_text(size=14), panel.border = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
well_temp_centroid
ggsave(well_temp_centroid, file="Graphs/multidirectional_temp_well_new.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
temp_radians <- temp_centroid_well2[,c("Common_name", "radians")]
colnames(temp_radians)[2] <- "temp_radians"
range_radians <- range_centroid_well2[,c("Common_name", "radians")]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(temp_radians, range_radians, by="Common_name", all.x=TRUE)
circ.cor(radians$temp_radians, radians$range_radians, test=TRUE)
## SIGNIFICANT (p<0.001) - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.5

# 2b: Mean total precipitation
#### Direction and distance of shift based on mean annual temperature
## 2.5th and 97.5th percentiles for each species in TP1
detach(package:Rmisc)
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_well %>% group_by(Common_name) %>% filter(Time_period=="1975-1991") %>%
  summarise(precip_2.5_perc = quantile(total_precip, probs = 0.025),
            precip_97.5_perc = quantile(total_precip, probs = 0.975))
nmrsdata_well <- merge(nmrsdata_well, tp1_species_percentiles, by="Common_name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_well %>% group_by(Common_name, Time_period) %>%
  filter(total_precip >= precip_2.5_perc & total_precip <= precip_97.5_perc)

## calculate centroids of each
precip_centroid_well <- filtered %>% group_by(Common_name, Time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
precip_centroid_well2 <- precip_centroid_well %>%
  gather(key, value, -Common_name, -Time_period) %>%
  unite(col, key, Time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2
sp_all <- data.frame(unique(nmrsdata_well$Common_name))
sp_sub <- data.frame(unique(precip_centroid_well2$Common_name))
sp_list <- sp_all %>% 
  filter(!sp_all$unique.nmrsdata_well.Common_name. %in% sp_sub$unique.precip_centroid_well2.Common_name.)

## remove NAs
precip_centroid_well2 <- na.omit(precip_centroid_well2) ## 68 species (6 removed)

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
precip_centroid_well2$distance <- sqrt((precip_centroid_well2$`mean_northing_1975-1991`-precip_centroid_well2$`mean_northing_2012-2016`)^2 + 
                                        (precip_centroid_well2$`mean_easting_1975-1991`-precip_centroid_well2$`mean_easting_2012-2016`)^2) ## pythagoras equation
precip_centroid_well2$distance <- precip_centroid_well2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
precip_centroid_well2$bearing <- bearing(precip_centroid_well2[,c(6,4)], precip_centroid_well2[,c(7,5)])
precip_centroid_well2$direction <- (precip_centroid_well2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
precip_centroid_well2$radians <- NISTdegTOradian(precip_centroid_well2$direction)
x <- deg(circ.mean(precip_centroid_well2$radians)) # -36, add full circle
(x+360)%%360 ## 340.04
ci <- vm.bootstrap.ci(precip_centroid_well2$radians) 
deg(c(-1.6, 0.38)) ## upper=22
(-91.67325+360)%%360 ## lower=268
# median & quartiles for distances
median(precip_centroid_well2$distance) ## 33.84 km 
quantile(precip_centroid_well2$distance) ## 25th = 16.13km, 75th = 53.93km

## plot result
ybreaks <- c(0,50,150,250)
well_precip_centroid <- ggplot(data=precip_centroid_well2,
                               aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1), lwd=0.6) +
  #geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F) +
  theme_bw() +
  theme(text = element_text(size=14), panel.border = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
well_precip_centroid
ggsave(well_precip_centroid, file="Graphs/multidirectional_total_precip_well_new.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
precip_radians <- precip_centroid_well2[,c("Common_name", "radians")]
colnames(precip_radians)[2] <- "precip_radians"
range_radians <- range_centroid_well2[,c("Common_name", "radians")]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(precip_radians, range_radians, by="Common_name", all.x=TRUE)
circ.cor(radians$precip_radians, radians$range_radians, test=TRUE)
## significant - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.74

## put all circular plots together
library(ggpubr)
library(gridExtra)
library(grid)
circular_well <- ggarrange(well_centroid, well_temp_centroid, well_precip_centroid,
                          labels = c("(a)", "(b)", "(c)"),font.label = list(size = 12),
                          ncol=3, nrow=1)
circular_well
ggsave(circular_well, file="Outputs/Graphs/circlar_range_temp_precip_well_new.png", height=8, width=10)

## Put datasets together
range_centroid_well2 <- range_centroid_well2[,c("Common_name", "distance", "direction")]
colnames(range_centroid_well2) <- c("Common_name", "range_distance", "range_direction")
temp_centroid_well2 <- temp_centroid_well2[,c("Common_name", "distance", "direction")]
colnames(temp_centroid_well2) <- c("Common_name", "temp_distance", "temp_direction")
precip_centroid_well2 <- precip_centroid_well2[,c("Common_name", "distance", "direction")]
colnames(precip_centroid_well2) <- c("Common_name", "precip_distance", "precip_direction")

well_centroid_shifts <- merge(range_centroid_well2, temp_centroid_well2, by="Common_name", all=TRUE)
well_centroid_shifts <- merge(well_centroid_shifts, precip_centroid_well2, by="Common_name", all=TRUE)
write.csv(well_centroid_shifts, file="Outputs/Results/NMRS_well_direction_distance.csv", row.names=FALSE)


#########################################################################################################
#########################################################################################################

rm(list = ls())
options(scipen=999)

## NMRS data for cool-adapted moths with elevation and temperature data at a 10km scale
nmrsdata <- readRDS("Data/NMRS/NMRS_cool_moths_final.rds") ## NMRS data for all hectads and all years with elevation
## hectad recording levels
hec_records <- read.csv("Data/NMRS/Hectad_recording_levels_1975_1991_2012_2016.csv", header=TRUE)

nmrsdata <- nmrsdata[!nmrsdata$Common_name=="Lunar Thorn",]

###################################
## Heavily recorded hectads

## filter to recorded hectads only (i.e. those recorded once in both time periods)
heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
heavy_hecs <- heavy_hecs[which(heavy_hecs$n_row > 1), ]
heavy_hecs$n_row <- NULL ## 1084 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_heavy <- nmrsdata[which(nmrsdata$Hectad %in% heavy_hecs$HECTAD), ]
length(unique(nmrsdata_heavy$Hectad)) ## 592 hectads
length(unique(nmrsdata_heavy$Common_name)) ## 72 species

#### 1. Multidimensional RANGE shift 
# Distance and direction of centroid shift between TP1 and TP2
# take mean easting and northing and lat and lon across each species occupied hectads for each time period = range centroid
# following Gillings et al 2015
range_centroid_heavy <- nmrsdata_heavy %>% group_by(Common_name, Time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
range_centroid_heavy2 <- range_centroid_heavy %>%
  gather(key, value, -Common_name, -Time_period) %>%
  unite(col, key, Time_period) %>%
  spread(col, value)
## remove NAs - some species are only found in one TP - Grey and Scotch Burnet removed
range_centroid_heavy2 <- na.omit(range_centroid_heavy2) ## 71 species
# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
range_centroid_heavy2$distance <- sqrt((range_centroid_heavy2$`mean_northing_1975-1991`-range_centroid_heavy2$`mean_northing_2012-2016`)^2 + 
                                        (range_centroid_heavy2$`mean_easting_1975-1991`-range_centroid_heavy2$`mean_easting_2012-2016`)^2) ## pythagoras equation
range_centroid_heavy2$distance <- range_centroid_heavy2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
range_centroid_heavy2$bearing <- bearing(range_centroid_heavy2[,c(6,4)], range_centroid_heavy2[,c(7,5)])
range_centroid_heavy2$direction <- (range_centroid_heavy2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360

## plot result
ybreaks <- c(0,100,200,300)
heavy_centroid <- ggplot(data=range_centroid_heavy2,
                        aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1), lwd=0.6) +
  #geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F) +
  theme_bw() +
  theme(text = element_text(size=14), panel.border = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
heavy_centroid
ggsave(heavy_centroid, file="Graphs/multidirectional_shifts_heavy.png")
## Welsh Clearwing south
## Sword-grass north

# use circ.mean from the circstats package to find a mean direction & confidence limits across all species
library(NISTunits)
library(CircStats)
range_centroid_heavy2$radians <- NISTdegTOradian(range_centroid_heavy2$direction)
x <- deg(circ.mean(range_centroid_heavy2$radians)) # -36, add full circle
(x+360)%%360 ## 331
ci <- vm.bootstrap.ci(range_centroid_heavy2$radians) 
deg(c(-3.63,0.76)) ## lower=152 upper=44
((-207.98368)+360)%%360 

# median & quartiles for distances
median(range_centroid_heavy2$distance) ## 41km 
quantile(range_centroid_heavy2$distance) ## 25th = 18km, 75th = 61km

## Rayleigh's test
## if significant, the distribution of directions differ significantly from an even distribution
## if R = 0 (the data are completely spread around the circle)
## if R = 1 (the data are completely concentrated on one point)
r.test(range_centroid_heavy2$radians)
## p = 0.37
## r = 0.12
## no significant bias in any direction


#### 2. Multidimensional CLIMATE shift
## What direction should species move in to remain in the same climatic envelope? 
# Then match this with the direction and distance species have actually shifted
# Do they correlate? 

##############
# 2a: Mean annual temperature
##############

#### Direction and distance of shift based on mean annual temperature
## 2.5th and 97.5th percentiles for each species in TP1
detach(package:Rmisc)
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_heavy %>% group_by(Common_name) %>% filter(Time_period=="1975-1991") %>%
  dplyr::summarise(temp_2.5_perc = quantile(temperature, probs = 0.025),
                   temp_97.5_perc = quantile(temperature, probs = 0.975))
nmrsdata_heavy <- merge(nmrsdata_heavy, tp1_species_percentiles, by="Common_name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_heavy %>% group_by(Common_name, Time_period) %>%
  filter(temperature >= temp_2.5_perc & temperature <= temp_97.5_perc)

## calculate centroids of each
temp_centroid_heavy <- filtered %>% group_by(Common_name, Time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing)) # 74 species

## change orientation of data
temp_centroid_heavy2 <- temp_centroid_heavy %>%
  gather(key, value, -Common_name, -Time_period) %>%
  unite(col, key, Time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2
sp_all <- data.frame(unique(nmrsdata_heavy$Common_name))
sp_sub <- data.frame(unique(temp_centroid_heavy2$Common_name))
sp_list <- sp_all %>% 
  filter(!sp_all$unique.nmrsdata_heavy.Common_name. %in% sp_sub$unique.temp_centroid_heavy2.Common_name.)

## remove NAs
temp_centroid_heavy2 <- na.omit(temp_centroid_heavy2) ## 67 species (7 species removed)

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
temp_centroid_heavy2$distance <- sqrt((temp_centroid_heavy2$`mean_northing_1975-1991`-temp_centroid_heavy2$`mean_northing_2012-2016`)^2 + 
                                       (temp_centroid_heavy2$`mean_easting_1975-1991`-temp_centroid_heavy2$`mean_easting_2012-2016`)^2) ## pythagoras equation
temp_centroid_heavy2$distance <- temp_centroid_heavy2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
temp_centroid_heavy2$bearing <- bearing(temp_centroid_heavy2[,c(6,4)], temp_centroid_heavy2[,c(7,5)])
temp_centroid_heavy2$direction <- (temp_centroid_heavy2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
temp_centroid_heavy2$radians <- NISTdegTOradian(temp_centroid_heavy2$direction)
x <- deg(circ.mean(temp_centroid_heavy2$radians)) # -11, add full circle
(x+360)%%360 ## 347
ci <- vm.bootstrap.ci(temp_centroid_heavy2$radians) 
deg(c(-0.4, -0.06)) 
(-22.918312+360)%%360 ## low = 337
(-3.437747+360)%%360 ## high = 357

# median & quartiles for distances
median(temp_centroid_heavy2$distance) ## 72km
quantile(temp_centroid_heavy2$distance) ## 25th = 41km, 75th = 107km

## plot result
ybreaks <- c(0,50,150,250)
heavy_temp_centroid <- ggplot(data=temp_centroid_heavy2,
                             aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1), lwd=0.6) +
  #geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F) +
  theme_bw() +
  theme(text = element_text(size=14), panel.border = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
heavy_temp_centroid
ggsave(heavy_temp_centroid, file="Graphs/multidirectional_temp_heavy.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
temp_radians <- temp_centroid_heavy2[,c("Common_name", "radians")]
colnames(temp_radians)[2] <- "temp_radians"
range_radians <- range_centroid_heavy2[,c("Common_name", "radians")]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(temp_radians, range_radians, by="Common_name", all.x=TRUE)
circ.cor(radians$temp_radians, radians$range_radians, test=TRUE)
## SIGNIFICANT (p=0.014) - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.36
## significant but not very correlated

# 2b: Mean total precipitation
#### Direction and distance of shift based on mean annual temperature
## 2.5th and 97.5th percentiles for each species in TP1
detach(package:Rmisc)
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_heavy %>% group_by(Common_name) %>% filter(Time_period=="1975-1991") %>%
  summarise(precip_2.5_perc = quantile(total_precip, probs = 0.025),
            precip_97.5_perc = quantile(total_precip, probs = 0.975))
nmrsdata_heavy <- merge(nmrsdata_heavy, tp1_species_percentiles, by="Common_name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_heavy %>% group_by(Common_name, Time_period) %>%
  filter(total_precip >= precip_2.5_perc & total_precip <= precip_97.5_perc)

## calculate centroids of each
precip_centroid_heavy <- filtered %>% group_by(Common_name, Time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
precip_centroid_heavy2 <- precip_centroid_heavy %>%
  gather(key, value, -Common_name, -Time_period) %>%
  unite(col, key, Time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2
sp_all <- data.frame(unique(nmrsdata_heavy$Common_name))
sp_sub <- data.frame(unique(precip_centroid_heavy2$Common_name))
sp_list <- sp_all %>% 
  filter(!sp_all$unique.nmrsdata_heavy.Common_name. %in% sp_sub$unique.precip_centroid_heavy2.Common_name.)

## remove NAs
precip_centroid_heavy2 <- na.omit(precip_centroid_heavy2) ## 60 species

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
precip_centroid_heavy2$distance <- sqrt((precip_centroid_heavy2$`mean_northing_1975-1991`-precip_centroid_heavy2$`mean_northing_2012-2016`)^2 + 
                                         (precip_centroid_heavy2$`mean_easting_1975-1991`-precip_centroid_heavy2$`mean_easting_2012-2016`)^2) ## pythagoras equation
precip_centroid_heavy2$distance <- precip_centroid_heavy2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
precip_centroid_heavy2$bearing <- bearing(precip_centroid_heavy2[,c(6,4)], precip_centroid_heavy2[,c(7,5)])
precip_centroid_heavy2$direction <- (precip_centroid_heavy2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
precip_centroid_heavy2$radians <- NISTdegTOradian(precip_centroid_heavy2$direction)
x <- deg(circ.mean(precip_centroid_heavy2$radians)) # -36, add full circle
(x+360)%%360 ## 32
ci <- vm.bootstrap.ci(precip_centroid_heavy2$radians) 
deg(c( -0.36, 2.6)) ## upper=149
(-20.62648+360)%%360 ## lower=339
# median & quartiles for distances
median(precip_centroid_heavy2$distance) ## 41 km 
quantile(precip_centroid_heavy2$distance) ## 25th = 20km, 75th = 61km

## plot result
ybreaks <- c(0,50,100,150,200,250,300)
heavy_precip_centroid <- ggplot(data=precip_centroid_heavy2,
                               aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1), lwd=0.6) +
  #geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F) +
  theme_bw() +
  theme(text = element_text(size=14), panel.border = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
heavy_precip_centroid
ggsave(heavy_precip_centroid, file="Graphs/multidirectional_total_precip_heavy.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
precip_radians <- precip_centroid_heavy2[,c("Common_name", "radians")]
colnames(precip_radians)[2] <- "precip_radians"
range_radians <- range_centroid_heavy2[,c("Common_name", "radians")]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(precip_radians, range_radians, by="Common_name", all.x=TRUE)
circ.cor(radians$precip_radians, radians$range_radians, test=TRUE)
## significant - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.3, p=0.014

## put all circular plots together
library(ggpubr)
library(gridExtra)
library(grid)
circular_heavy <- ggarrange(heavy_centroid, heavy_temp_centroid, heavy_precip_centroid,
                           labels = c("(a)", "(b)", "(c)"),font.label = list(size = 12),
                           ncol=3, nrow=1)
circular_heavy
ggsave(circular_heavy, file="Outputs/Graphs/circlar_range_temp_precip_heavy.png", height=8, width=10)

## Put datasets together
range_centroid_heavy2 <- range_centroid_heavy2[,c("Common_name", "distance", "direction")]
colnames(range_centroid_heavy2) <- c("Common_name", "range_distance", "range_direction")
temp_centroid_heavy2 <- temp_centroid_heavy2[,c("Common_name", "distance", "direction")]
colnames(temp_centroid_heavy2) <- c("Common_name", "temp_distance", "temp_direction")
precip_centroid_heavy2 <- precip_centroid_heavy2[,c("Common_name", "distance", "direction")]
colnames(precip_centroid_heavy2) <- c("Common_name", "precip_distance", "precip_direction")

heavy_centroid_shifts <- merge(range_centroid_heavy2, temp_centroid_heavy2, by="Common_name", all=TRUE)
heavy_centroid_shifts <- merge(heavy_centroid_shifts, precip_centroid_heavy2, by="Common_name", all=TRUE)

write.csv(heavy_centroid_shifts, file="Outputs/Results/NMRS_heavy_direction_distance.csv", row.names=FALSE)
heavy_centroid_shifts$recording_effort <- "Heavily recorded"
well_centroid_shifts$recording_effort <- "Well recorded"

all_centroid_shifts <- rbind(well_centroid_shifts, heavy_centroid_shifts)
write.csv(all_centroid_shifts, file="Outputs/Results/NMRS_well_heavy_direction_distance.csv", row.names=FALSE)

nmrs_direction <- read.csv("Outputs/Results/NMRS_well_heavy_direction_distance.csv", header=TRUE)
# add in grey - can't be calculated for any centroid shift as no records in TP2 but still a cool-adapted species
library(tidyverse)
nmrs_direction <- nmrs_direction %>% add_row(Common_name = "Grey", recording_effort="Well recorded")
## populate so each species has a row for each recording effort level
nmrs_direction2 <- nmrs_direction %>%
  expand(Common_name, recording_effort)
nmrs_direction2 <- merge(nmrs_direction2, nmrs_direction, by=c("Common_name", "recording_effort"), all.x=TRUE)
sci_names <- unique(nmrsdata_well[,c("Common_name", "Scientific_name")]) ## add scientific names
nmrs_direction2 <- merge(nmrs_direction2, sci_names, by="Common_name", all.x=TRUE)
nmrs_direction2 <- nmrs_direction2[order(nmrs_direction2$Scientific_name),]
nmrs_direction2 <- nmrs_direction2 %>% group_by(Scientific_name) %>% arrange(desc(recording_effort), .by_group = TRUE)
nmrs_direction2 <- nmrs_direction2[,c(9,1:8)]
write.csv(nmrs_direction2, file="Outputs/Results/NMRS_well_heavy_direction_distance_final.csv", row.names=FALSE)

#
