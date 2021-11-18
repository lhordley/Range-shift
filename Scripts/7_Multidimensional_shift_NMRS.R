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

###############################################################
###################### RECORDED HECTADS #######################
###############################################################

## filter to recorded hectads only (i.e. those recorded once in both time periods)
hectads <- unique(hec_records$HECTAD) ## 1773
nmrsdata_rec <- nmrsdata[which(nmrsdata$Hectad %in% hectads), ]
length(unique(nmrsdata_rec$Hectad)) ## 1424 hectads
length(unique(nmrsdata_rec$Common_Name)) ## 72 species

#### 1. Multidimensional RANGE shift 
# Distance and direction of centroid shift between TP1 and TP2
# take mean easting and northing and lat and lon across each species occupied hectads for each time period = range centroid
# following Gillings et al 2015
range_centroid_rec <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
range_centroid_rec2 <- range_centroid_rec %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
range_centroid_rec2$distance <- sqrt((range_centroid_rec2$mean_northing_TP1-range_centroid_rec2$mean_northing_TP2)^2 + 
                                       (range_centroid_rec2$mean_easting_TP1-range_centroid_rec2$mean_easting_TP2)^2) ## pythagoras equation
range_centroid_rec2$distance <- range_centroid_rec2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
range_centroid_rec2$bearing <- bearing(range_centroid_rec2[,c(6,4)], range_centroid_rec2[,c(7,5)])
range_centroid_rec2$direction <- (range_centroid_rec2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360

## plot result
ybreaks <- c(0,50,100,150)
rec_centroid <- ggplot(data=range_centroid_rec2,
                       aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F,
            size = 3) +
  theme_bw() +
 theme( panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank())
rec_centroid
ggsave(rec_centroid, file="Graphs/multidirectional_shifts_rec.png")

# use circ.mean from the circstats package to find a mean direction & confidence limits across all species
library(NISTunits)
library(CircStats)
range_centroid_rec2$radians <- NISTdegTOradian(range_centroid_rec2$direction)
x <- deg(circ.mean(range_centroid_rec2$radians)) # -36, add full circle
(x+360)%%360 ## 323
# median & quartiles for distances
median(range_centroid_rec2$distance) ## 40.2km 
quantile(range_centroid_rec2$distance) ## 25th = 22.9km, 75th = 70.9km

## Rayleigh's test
## if significant, the distribution of directions differ significantly from an even distribution
## if R = 0 (the data are completely spread around the circle)
## if R = 1 (the data are completely concentrated on one point)
r.test(range_centroid_rec2$radians)
## p < 0.001
## r = 0.42
## low R = variable directions, but significant means it has significant bias to a certain direction (the mean)

## Are species that are moving southwards moving uphill?
# Test by calculating mean elevation change between time periods
mean_elev_rec <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_elev=mean(elevation10x10km))
## change orientation
mean_elev_rec2 <- mean_elev_rec %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
mean_elev_rec2$elev_shift <- mean_elev_rec2$mean_elev_TP2 - mean_elev_rec2$mean_elev_TP1 ## TP2 - TP1
hist(mean_elev_rec2$elev_shift) ## looks good

mean(mean_elev_rec2$elev_shift) ## 38.23 m uphill
sd(mean_elev_rec2$elev_shift) ## 111.27m 

# Shapiro-Wilk normality test for the differences
shapiro.test(mean_elev_rec2$elev_shift) # => p-value = 0.02 - distribution is NOT normal
qqnorm(mean_elev_rec2$elev_shift)
qqline(mean_elev_rec2$elev_shift) ## not good

## Wilcox test to test whether elevation shifts differ significantly from 0
wilcox.test(mean_elev_rec2$elev_shift) ## significant

## wilcoxon signed-rank test - does mean elevation differ between time periods?
elev_rec <- wilcox.test(mean_elev ~ time_period, data = mean_elev_rec, paired = TRUE)
elev_rec ## significant

## plot result
elev_rec_plot <- ggpaired(mean_elev_rec, x = "time_period", y = "mean_elev",
                           color = "time_period", line.color = "gray", line.size = 0.4,
                           palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Mean elevation")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
elev_rec_plot ## higher mean elevation in TP2 compared to TP1
ggsave(elev_rec_plot, file="Graphs/Mean_elevation_shift_rec_hecs.png")

length(which(mean_elev_rec2$elev_shift>0)) ## 49 species moving uphill
length(which(mean_elev_rec2$elev_shift<0)) ## 23 species moving downhill
## categorise species as either moving uphill or downhill
mean_elev_rec2$shift_direction <- ifelse(mean_elev_rec2$elev_shift>0, "Uphill", "Downhill")

## add in elevation results to direction and distance shifts
centroid_elev_rec <- merge(mean_elev_rec2, range_centroid_rec2, by="Common_Name")

## plot result coloured by elevational shift
rec_centroid2 <- ggplot(data=centroid_elev_rec,
                       aes(x=direction, y=distance, colour=shift_direction)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  ggtitle("Recorded hectads") +
  theme_bw()
rec_centroid2
ggsave(rec_centroid2, file="Graphs/multidirectional_elev_shifts_rec.png")


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
tp1_species_percentiles <- nmrsdata_rec %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
  summarise(temp_2.5_perc = quantile(temperature, probs = 0.025),
            temp_97.5_perc = quantile(temperature, probs = 0.975))
nmrsdata_rec <- merge(nmrsdata_rec, tp1_species_percentiles, by="Common_Name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>%
  filter(temperature >= temp_2.5_perc & temperature <= temp_97.5_perc)

## calculate centroids of each
temp_centroid_rec <- filtered %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
temp_centroid_rec2 <- temp_centroid_rec %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2

## remove NAs
temp_centroid_rec2 <- na.omit(temp_centroid_rec2) ## 65 species (7 removed)

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
temp_centroid_rec2$distance <- sqrt((temp_centroid_rec2$mean_northing_TP1-temp_centroid_rec2$mean_northing_TP2)^2 + 
                                      (temp_centroid_rec2$mean_easting_TP1-temp_centroid_rec2$mean_easting_TP2)^2) ## pythagoras equation
temp_centroid_rec2$distance <- temp_centroid_rec2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
temp_centroid_rec2$bearing <- bearing(temp_centroid_rec2[,c(6,4)], temp_centroid_rec2[,c(7,5)])
temp_centroid_rec2$direction <- (temp_centroid_rec2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
temp_centroid_rec2$radians <- NISTdegTOradian(temp_centroid_rec2$direction)
x <- deg(circ.mean(temp_centroid_rec2$radians)) # -36, add full circle
(x+360)%%360 ## 337
# median & quartiles for distances
median(temp_centroid_rec2$distance) ## 60.06km 
quantile(temp_centroid_rec2$distance) ## 25th = 30.5km, 75th = 93.4km

## plot result
ybreaks <- c(0,50,100,150,200,250,300)
rec_temp_centroid <- ggplot(data=temp_centroid_rec2,
                       aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) +
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F,
            size = 3) +
  theme_bw() +
  theme( panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
         axis.title.x = element_blank(), axis.title.y = element_blank())
rec_temp_centroid
ggsave(rec_temp_centroid, file="Graphs/multidirectional_temp_rec.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
temp_radians <- temp_centroid_rec2[,c(1,13)]
colnames(temp_radians)[2] <- "temp_radians"
range_radians <- range_centroid_rec2[,c(1,13)]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(temp_radians, range_radians, by="Common_Name", all.x=TRUE)
circ.cor(radians$temp_radians, radians$range_radians, test=TRUE)
## significant - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.749

###################
# 2b: Mean summer temperature
###################

#### Direction and distance of shift based on mean annual temperature
## 2.5th and 97.5th percentiles for each species in TP1
detach(package:Rmisc)
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_rec %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
  summarise(summer_temp_2.5_perc = quantile(summer_temperature, probs = 0.025),
            summer_temp_97.5_perc = quantile(summer_temperature, probs = 0.975))
nmrsdata_rec <- merge(nmrsdata_rec, tp1_species_percentiles, by="Common_Name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>%
  filter(summer_temperature >= summer_temp_2.5_perc & summer_temperature <= summer_temp_97.5_perc)

## calculate centroids of each
summer_temp_centroid_rec <- filtered %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
summer_temp_centroid_rec2 <- summer_temp_centroid_rec %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2

## remove NAs
summer_temp_centroid_rec2 <- na.omit(summer_temp_centroid_rec2) ## 66 species (7 removed)

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
summer_temp_centroid_rec2$distance <- sqrt((summer_temp_centroid_rec2$mean_northing_TP1-summer_temp_centroid_rec2$mean_northing_TP2)^2 + 
                                      (summer_temp_centroid_rec2$mean_easting_TP1-summer_temp_centroid_rec2$mean_easting_TP2)^2) ## pythagoras equation
summer_temp_centroid_rec2$distance <- summer_temp_centroid_rec2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
summer_temp_centroid_rec2$bearing <- bearing(summer_temp_centroid_rec2[,c(6,4)], summer_temp_centroid_rec2[,c(7,5)])
summer_temp_centroid_rec2$direction <- (summer_temp_centroid_rec2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
summer_temp_centroid_rec2$radians <- NISTdegTOradian(summer_temp_centroid_rec2$direction)
x <- deg(circ.mean(summer_temp_centroid_rec2$radians)) # -36, add full circle
(x+360)%%360 ## 327
# median & quartiles for distances
median(summer_temp_centroid_rec2$distance) ## 46.3km 
quantile(summer_temp_centroid_rec2$distance) ## 25th = 27.18km, 75th = 72.8km

## plot result
ybreaks <- c(0,50,100,150,200,250,300)
rec_summer_temp_centroid <- ggplot(data=summer_temp_centroid_rec2,
                            aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F,
            size = 3) +
  theme_bw() +
  theme( panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
         axis.title.x = element_blank(), axis.title.y = element_blank())
rec_summer_temp_centroid
ggsave(rec_summer_temp_centroid, file="Graphs/multidirectional_summer_temp_rec.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
temp_radians <- summer_temp_centroid_rec2[,c(1,13)]
colnames(temp_radians)[2] <- "temp_radians"
range_radians <- range_centroid_rec2[,c(1,13)]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(temp_radians, range_radians, by="Common_Name", all.x=TRUE)
circ.cor(radians$temp_radians, radians$range_radians, test=TRUE)
## significant - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.752
## very similar to mean annual temperature

# 2c: Mean total precipitation
#### Direction and distance of shift based on mean annual temperature
## 2.5th and 97.5th percentiles for each species in TP1
detach(package:Rmisc)
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_rec %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
  summarise(precip_2.5_perc = quantile(total_precip, probs = 0.025),
            precip_97.5_perc = quantile(total_precip, probs = 0.975))
nmrsdata_rec <- merge(nmrsdata_rec, tp1_species_percentiles, by="Common_Name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>%
  filter(total_precip >= precip_2.5_perc & total_precip <= precip_97.5_perc)

## calculate centroids of each
precip_centroid_rec <- filtered %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
precip_centroid_rec2 <- precip_centroid_rec %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2

## remove NAs
precip_centroid_rec2 <- na.omit(precip_centroid_rec2) ## 67 species (5 removed)

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
precip_centroid_rec2$distance <- sqrt((precip_centroid_rec2$mean_northing_TP1-precip_centroid_rec2$mean_northing_TP2)^2 + 
                                      (precip_centroid_rec2$mean_easting_TP1-precip_centroid_rec2$mean_easting_TP2)^2) ## pythagoras equation
precip_centroid_rec2$distance <- precip_centroid_rec2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
precip_centroid_rec2$bearing <- bearing(precip_centroid_rec2[,c(6,4)], precip_centroid_rec2[,c(7,5)])
precip_centroid_rec2$direction <- (precip_centroid_rec2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
precip_centroid_rec2$radians <- NISTdegTOradian(precip_centroid_rec2$direction)
x <- deg(circ.mean(precip_centroid_rec2$radians)) # -36, add full circle
(x+360)%%360 ## 323
x.bs <- vm.bootstrap.ci(precip_centroid_rec2$radians)
(0.35+360)%%360 ## 323

library(Directional)
circ.summary(precip_centroid_rec2$radians, rads=TRUE)

r.test(precip_centroid_rec2$radians)

# median & quartiles for distances
median(precip_centroid_rec2$distance) ## 35.26km 
quantile(precip_centroid_rec2$distance) ## 25th = 20.12km, 75th = 57.12km

## plot result
ybreaks <- c(0,50,100,150,200,250,300)
rec_precip_centroid <- ggplot(data=precip_centroid_rec2,
                                   aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F,
            size = 3) +
  theme_bw() +
  theme( panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
         axis.title.x = element_blank(), axis.title.y = element_blank())
rec_precip_centroid
ggsave(rec_precip_centroid, file="Graphs/multidirectional_total_precip_rec.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
precip_radians <- precip_centroid_rec2[,c(1,13)]
colnames(precip_radians)[2] <- "precip_radians"
range_radians <- range_centroid_rec2[,c(1,13)]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(precip_radians, range_radians, by="Common_Name", all.x=TRUE)
circ.cor(radians$precip_radians, radians$range_radians, test=TRUE)
## significant - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.88
## very similar to mean annual temperature

## correlate the distances too
precip_dist <- precip_centroid_rec2[,c(1,10)]
colnames(precip_dist)[2] <- "precip_dist"
range_dist <- range_centroid_rec2[,c(1,10)]
colnames(range_dist)[2] <- "range_dist"
distance <- merge(precip_dist, range_dist, by="Common_Name", all.x=TRUE)

plot(distance$precip_dist, distance$range_dist)
cor.test(distance$precip_dist, distance$range_dist, method="spearman")
## 0.86, p<0.001


## put all circular plots together
library(ggpubr)
circular_rec <- ggarrange(rec_centroid, rec_temp_centroid, rec_summer_temp_centroid, rec_precip_centroid,
                          labels = c("(a)", "(b)", "(c)", "(d)"),
                          ncol = 2, nrow = 2)
ggsave(circular_rec, file="Graphs/circlar_range_temp_summer_precip_rec.png", height=12, width=10)

#### 3. Define colonised, persisted, and extirpated hectads
nmrsdata_rec$Recorded <- 1
nmrsdata_rec_expand <- nmrsdata_rec %>% expand(Common_Name, time_period, Hectad)
recorded <- nmrsdata_rec[,c("Common_Name", "time_period", "Hectad", "Recorded")]
nmrsdata_rec_expand <- merge(nmrsdata_rec_expand, recorded, by=c("Common_Name", "time_period", "Hectad"), all.x=TRUE)
nmrsdata_rec_expand[is.na(nmrsdata_rec_expand)] <- 0 ## change NAs to zero == species was NOT recorded at this hectad in this time period
## change to long format
nmrsdata_rec_expand2 <- nmrsdata_rec_expand %>%
  spread(time_period, Recorded) ## each species has 1424 rows = the number of recorded hectads
## first remove rows where TP1 AND TP2 == 0 (this a hectad where a species was never recorded)
nmrsdata_rec_expand2<-nmrsdata_rec_expand2[!(nmrsdata_rec_expand2$TP1==0 & nmrsdata_rec_expand2$TP2==0),]
nmrsdata_rec_expand2$Hectad_category <- case_when(
  nmrsdata_rec_expand2$TP1==0 & nmrsdata_rec_expand2$TP2==1 ~ "Colonisation",
  nmrsdata_rec_expand2$TP1==1 & nmrsdata_rec_expand2$TP2==0 ~ "Extirpation",
  TRUE ~ "Persistence"
  )

## plot results on a map
## add in lat/lon info
lat_lon <- unique(nmrsdata_rec[,c("Hectad", "lat", "lon","elevation10x10km","elevation10x10km_SD")])
nmrsdata_rec_expand2 <- merge(nmrsdata_rec_expand2, lat_lon, by="Hectad", all.x=TRUE)
worldmap = map_data('world')
ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  theme_void() +
  geom_point(data =nmrsdata_rec_expand2,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat), colour=Hectad_category), size=2) +
  theme(title = element_text(size = 12))

## find proportion of species at each extirpated hectad - where are extirpations occurring? 
extirpated_hecs1 <- nmrsdata_rec_expand2 %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
                   dplyr::summarise(tot_sp = n())
extirpated_hecs2 <- nmrsdata_rec_expand2 %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
  filter(Hectad_category=="Extirpation") %>%
  dplyr::summarise(extir_sp = n())
extirpated_hecs <- merge(extirpated_hecs1, extirpated_hecs2, by=c("Hectad", "lat", "lon", "elevation10x10km", "elevation10x10km_SD"))
extirpated_hecs$extir_prop <- extirpated_hecs$extir_sp/extirpated_hecs$tot_sp
## very long way of doing this - gave up on finding a quicker way!
                  
## heatmap
extir_hecs <- ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  theme_void() +
  geom_point(data =extirpated_hecs,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat), colour=extir_prop), size=1) +
  scale_color_viridis_c(name="Proportion of species") + 
  theme(title = element_text(size = 12))
## extirpated hectads with the most species lost tend to be in Northern Scotland/Pennines
extir_hecs
ggsave(extir_hecs, file="Maps/Rec_extirpated_hecs2.png")


#################################################################
#### How does the climate and elevation of extirpated hectads compared to all recorded hectads?

## proportion plots that Andy suggested
rec_hecs <- nmrsdata_rec[,c("Hectad","elevation10x10km","elevation10x10km_SD")] ## all rec hecs: 1424
rec_hecs <- rec_hecs %>% distinct(Hectad, .keep_all = TRUE) ## different way of unique
## subset 'well' extirpated hecs which have lost at least 25% of species
well_extir_hecs <- extirpated_hecs[extirpated_hecs$extir_prop>=0.25,] ## 814 hectads

## add in climate data to each data frame
temp <- read.csv("Data/NMRS/All_NMRS_hectads_annual_temperature.csv", header=TRUE)
summer_temp <- read.csv("Data/NMRS/All_NMRS_hectads_summer_temperature.csv", header=TRUE)
precip <- read.csv("Data/NMRS/All_NMRS_hectads_total_precipitation.csv", header=TRUE)
## put these together
df_list <- list(temp, summer_temp, precip)
climate <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)
## only want TP2 climate data
climate_tp2 <- climate[climate$time_period=="TP2",]

## merge in with all recorded hectads and then all well extirpated hectads
rec_hecs <- merge(rec_hecs, climate_tp2, by="Hectad", all.x=)
well_extir_hecs <- merge(well_extir_hecs, climate_tp2, by="Hectad", all.x=)

#### Annual temperature
# split temperature into 10 bins
rec_hecs$temperature_cat <- cut(rec_hecs$temperature, breaks = 10, labels=FALSE)

rec_hecs_temp_prop <- rec_hecs %>% group_by(temperature_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(rec_hecs_temp_prop$freq ~ rec_hecs_temp_prop$temperature_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(temperature_cat = cut(rec_hecs$temperature, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(temperature_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)
min_max$cat <- 1:10

medians <- rec_hecs %>% group_by(temperature_cat) %>% dplyr::summarise(median=median(temperature))

extirpated_hecs2 <- well_extir_hecs %>% mutate(temperature_cat = case_when(
  between(temperature, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(temperature, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(temperature, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(temperature, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(temperature, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(temperature, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(temperature, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(temperature, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(temperature, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(temperature_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$temperature_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("temperature_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
rec_hecs_temp_prop <- rec_hecs_temp_prop[,c("temperature_cat", "freq")]
rec_hecs_temp_prop$hec_cat <- "All hectads"

temp_prop <- rbind(extirpated_hecs_prop, rec_hecs_temp_prop)
temp_prop <- merge(temp_prop, medians, by="temperature_cat")
## plot line graph
rec_extir_temp <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Mean annual temperature", y="Proportion of hectads") +
  scale_x_continuous(limits=c(4,12)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
rec_extir_temp
ggsave(rec_extir_temp, file="Graphs/Rec_hecs_extirpated_temperature.png")


#### Summer temperature
# split temperature into 10 bins
rec_hecs$summer_temperature_cat <- cut(rec_hecs$summer_temperature, breaks = 10, labels=FALSE)

rec_hecs_temp_prop <- rec_hecs %>% group_by(summer_temperature_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(rec_hecs_temp_prop$freq ~ rec_hecs_temp_prop$summer_temperature_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(summer_temperature_cat = cut(rec_hecs$summer_temperature, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(summer_temperature_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)

medians <- rec_hecs %>% group_by(summer_temperature_cat) %>% dplyr::summarise(median=median(summer_temperature))

extirpated_hecs2 <- well_extir_hecs %>% mutate(summer_temperature_cat = case_when(
  between(summer_temperature, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(summer_temperature, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(summer_temperature, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(summer_temperature, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(summer_temperature, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(summer_temperature, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(summer_temperature, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(summer_temperature, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(summer_temperature, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(summer_temperature_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$summer_temperature_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("summer_temperature_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
rec_hecs_temp_prop <- rec_hecs_temp_prop[,c("summer_temperature_cat", "freq")]
rec_hecs_temp_prop$hec_cat <- "All hectads"

temp_prop <- rbind(extirpated_hecs_prop, rec_hecs_temp_prop)
temp_prop <- merge(temp_prop, medians, by="summer_temperature_cat")
## plot line graph
rec_extir_summer_temp <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Mean annual summer temperature", y="Proportion of hectads") +
  scale_x_continuous(limits=c(10,18)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
rec_extir_summer_temp
ggsave(rec_extir_summer_temp, file="Graphs/Rec_hecs_extirpated_summer_temperature.png")



#### Total precipitation
# split precipitatino into 10 bins
rec_hecs$precip_cat <- cut(rec_hecs$total_precip, breaks = 10, labels=FALSE)

rec_hecs_precip_prop <- rec_hecs %>% group_by(precip_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(rec_hecs_precip_prop$freq ~ rec_hecs_precip_prop$precip_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(precip_cat = cut(rec_hecs$total_precip, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(precip_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)

medians <- rec_hecs %>% group_by(precip_cat) %>% dplyr::summarise(median=median(total_precip))

extirpated_hecs2 <- well_extir_hecs %>% mutate(precip_cat = case_when(
  between(total_precip, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(total_precip, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(total_precip, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(total_precip, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(total_precip, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(total_precip, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(total_precip, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(total_precip, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(total_precip, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(precip_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$precip_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("precip_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
rec_hecs_precip_prop <- rec_hecs_precip_prop[,c("precip_cat", "freq")]
rec_hecs_precip_prop$hec_cat <- "All hectads"

precip_prop <- rbind(extirpated_hecs_prop, rec_hecs_precip_prop)
precip_prop <- merge(precip_prop, medians, by="precip_cat")
## plot line graph
rec_extir_precip <- ggplot(precip_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Mean total precipitation (mm)", y="Proportion of hectads") +
  scale_x_continuous(breaks=seq(500,3500, by=500)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
rec_extir_precip
ggsave(rec_extir_precip, file="Graphs/Rec_hecs_extirpated_total_precipitation.png")


#### Elevation
# split elevation into 10 bins
rec_hecs$elev_cat <- cut(rec_hecs$elevation10x10km, breaks = 10, labels=FALSE)

rec_hecs_elev_prop <- rec_hecs %>% group_by(elev_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(rec_hecs_elev_prop$freq ~ rec_hecs_elev_prop$elev_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(elev_cat = cut(rec_hecs$elevation10x10km, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(elev_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)

medians <- rec_hecs %>% group_by(elev_cat) %>% dplyr::summarise(median=median(elevation10x10km))

extirpated_hecs2 <- well_extir_hecs %>% mutate(elev_cat = case_when(
  between(elevation10x10km, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(elevation10x10km, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(elevation10x10km, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(elevation10x10km, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(elevation10x10km, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(elevation10x10km, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(elevation10x10km, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(elevation10x10km, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(elevation10x10km, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(elev_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$elev_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("elev_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
rec_hecs_elev_prop <- rec_hecs_elev_prop[,c("elev_cat", "freq")]
rec_hecs_elev_prop$hec_cat <- "All hectads"

elev_prop <- rbind(extirpated_hecs_prop, rec_hecs_elev_prop)
elev_prop <- merge(elev_prop, medians, by="elev_cat")
## plot line graph
rec_extir_elev <- ggplot(elev_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Mean elevation (m)", y="Proportion of hectads") +
  scale_x_continuous(breaks=seq(0,2000, by=200)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
rec_extir_elev
ggsave(rec_extir_elev, file="Graphs/Rec_hecs_extirpated_elevation.png")


#### SD Elevation
# split elevation into 10 bins
rec_hecs$elev_cat <- cut(rec_hecs$elevation10x10km_SD, breaks = 10, labels=FALSE)

rec_hecs_elev_prop <- rec_hecs %>% group_by(elev_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(rec_hecs_elev_prop$freq ~ rec_hecs_elev_prop$elev_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(elev_cat = cut(rec_hecs$elevation10x10km_SD, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(elev_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)

medians <- rec_hecs %>% group_by(elev_cat) %>% dplyr::summarise(median=median(elevation10x10km_SD))

extirpated_hecs2 <- well_extir_hecs %>% mutate(elev_cat = case_when(
  between(elevation10x10km_SD, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(elevation10x10km_SD, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(elevation10x10km_SD, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(elevation10x10km_SD, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(elevation10x10km_SD, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(elevation10x10km_SD, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(elevation10x10km_SD, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(elevation10x10km_SD, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(elevation10x10km_SD, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(elev_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$elev_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("elev_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
rec_hecs_elev_prop <- rec_hecs_elev_prop[,c("elev_cat", "freq")]
rec_hecs_elev_prop$hec_cat <- "All hectads"

elev_prop <- rbind(extirpated_hecs_prop, rec_hecs_elev_prop)
elev_prop <- merge(elev_prop, medians, by="elev_cat")
## plot line graph
rec_extir_elev_sd <- ggplot(elev_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Standard deviation elevation (m)", y="Proportion of hectads") +
  scale_x_continuous(limits=c(0,400)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
rec_extir_elev_sd
ggsave(rec_extir_elev_sd, file="Graphs/Rec_hecs_extirpated_elevation_SD.png")


### Are hectads which have warmed more over time more likely to be extirpated?
## calculate difference in climate (tp2 - tp1)
temp_df <- climate[,c("Hectad", "time_period", "temperature")]
climate_diff <- temp_df %>%
  spread(time_period, temperature) %>% 
  mutate(temp_diff = TP2-TP1)

## merge in with all recorded hectads and then all well extirpated hectads
rec_hecs <- merge(rec_hecs, climate_diff, by="Hectad", all.x=TRUE)
well_extir_hecs <- merge(well_extir_hecs, climate_diff, by="Hectad", all.x=TRUE)

#### Annual temperature
# split temperature into 10 bins
rec_hecs$temp_diff_cat <- cut(rec_hecs$temp_diff, breaks = 10, labels=FALSE)

rec_hecs_temp_prop <- rec_hecs %>% group_by(temp_diff_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(rec_hecs_temp_prop$freq ~ rec_hecs_temp_prop$temp_diff_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(temp_diff_cat = cut(rec_hecs$temp_diff, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(temp_diff_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)

medians <- rec_hecs %>% group_by(temp_diff_cat) %>% dplyr::summarise(median=median(temp_diff))

extirpated_hecs2 <- well_extir_hecs %>% mutate(temp_diff_cat = case_when(
  between(temp_diff, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(temp_diff, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(temp_diff, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(temp_diff, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(temp_diff, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(temp_diff, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(temp_diff, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(temp_diff, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(temp_diff, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(temp_diff_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$temp_diff_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("temp_diff_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
rec_hecs_temp_prop <- rec_hecs_temp_prop[,c("temp_diff_cat", "freq")]
rec_hecs_temp_prop$hec_cat <- "All hectads"

temp_prop <- rbind(extirpated_hecs_prop, rec_hecs_temp_prop)
temp_prop <- merge(temp_prop, medians, by="temp_diff_cat")
## plot line graph
rec_extir_temp_diff <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Difference in temperature", y="Proportion of hectads") +
  scale_x_continuous(limits=c(0.25,1.15)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
rec_extir_temp_diff
ggsave(rec_extir_temp_diff, file="Graphs/Rec_hecs_extirpated_temperature_difference.png")




##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

rm(list = ls())
options(scipen=999)

## NMRS data for cool-adapted moths with elevation and temperature data at a 10km scale
nmrsdata <- readRDS("Data/NMRS/NMRS_cool_moths_final.rds") ## NMRS data for all hectads and all years with elevation
## hectad recording levels
hec_records <- read.csv("Data/NMRS/Hectad_recording_levels_1975_1991_2012_2016.csv", header=TRUE)


###################################
## well & heavily recorded hectads

## filter to recorded hectads only (i.e. those recorded once in both time periods)
well_heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded" | RECORDING.LEVEL=="Well recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
well_heavy_hecs <- well_heavy_hecs[which(well_heavy_hecs$n_row > 1), ]
well_heavy_hecs$n_row <- NULL ## 747 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_well <- nmrsdata[which(nmrsdata$Hectad %in% well_heavy_hecs$HECTAD), ]
length(unique(nmrsdata_well$Hectad)) ## 639 hectads
length(unique(nmrsdata_well$Common_Name)) ## 69 species

#### 1. Multidimensional RANGE shift 
# Distance and direction of centroid shift between TP1 and TP2
# take mean easting and northing and lat and lon across each species occupied hectads for each time period = range centroid
# following Gillings et al 2015
range_centroid_well <- nmrsdata_well %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
range_centroid_well2 <- range_centroid_well %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
## remove NAs - some species are only found in one TP 
range_centroid_well2 <- na.omit(range_centroid_well2) ## 67 species
# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
range_centroid_well2$distance <- sqrt((range_centroid_well2$mean_northing_TP1-range_centroid_well2$mean_northing_TP2)^2 + 
                                       (range_centroid_well2$mean_easting_TP1-range_centroid_well2$mean_easting_TP2)^2) ## pythagoras equation
range_centroid_well2$distance <- range_centroid_well2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
range_centroid_well2$bearing <- bearing(range_centroid_well2[,c(6,4)], range_centroid_well2[,c(7,5)])
range_centroid_well2$direction <- (range_centroid_well2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360

## plot result
ybreaks <- c(0,50,100,150,200)
well_centroid <- ggplot(data=range_centroid_well2,
                       aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) +
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F,
            size = 3) +
  theme_bw() +
  theme( panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
         axis.title.x = element_blank(), axis.title.y = element_blank())
well_centroid
ggsave(well_centroid, file="Graphs/multidirectional_shifts_well.png")

# use circ.mean from the circstats package to find a mean direction & confidence limits across all species
library(NISTunits)
library(CircStats)
range_centroid_well2$radians <- NISTdegTOradian(range_centroid_well2$direction)
x <- deg(circ.mean(range_centroid_well2$radians)) # -36, add full circle
(x+360)%%360 ## 230
# median & quartiles for distances
median(range_centroid_well2$distance) ## 37.2km 
quantile(range_centroid_well2$distance) ## 25th = 22.8km, 75th = 68km

## Rayleigh's test
## if significant, the distribution of directions differ significantly from an even distribution
## if R = 0 (the data are completely spread around the circle)
## if R = 1 (the data are completely concentrated on one point)
r.test(range_centroid_well2$radians)
## p < 0.001
## r = 0.36
## low R = variable directions, but significant means it has significant bias to a certain direction (the mean)

## Are species that are moving southwards moving uphill?
# Test by calculating mean elevation change between time periods
mean_elev_well <- nmrsdata_well %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_elev=mean(elevation10x10km))
## change orientation
mean_elev_well2 <- mean_elev_well %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
## remove NAs - some species are only found in one TP 
mean_elev_well2 <- na.omit(mean_elev_well2) ## 67 species

mean_elev_well2$elev_shift <- mean_elev_well2$mean_elev_TP2 - mean_elev_well2$mean_elev_TP1 ## TP2 - TP1
hist(mean_elev_well2$elev_shift) ## looks good

mean(mean_elev_well2$elev_shift) ## 29.4 m uphill
sd(mean_elev_well2$elev_shift) ## 142.5m 

# Shapiro-Wilk normality test for the differences
shapiro.test(mean_elev_well2$elev_shift) # => p-value = 0.02 - distribution is NOT normal
qqnorm(mean_elev_well2$elev_shift)
qqline(mean_elev_well2$elev_shift) ## not good

## Wilcox test to test whether elevation shifts differ significantly from 0
wilcox.test(mean_elev_well2$elev_shift) ## significant

## wilcoxon signed-rank test - does mean elevation differ between time periods?
elev_well <- wilcox.test(mean_elev ~ time_period, data = mean_elev_well, paired = TRUE)
elev_well ## non-significant

## plot result
elev_well_plot <- ggpaired(mean_elev_well, x = "time_period", y = "mean_elev",
                          color = "time_period", line.color = "gray", line.size = 0.4,
                          palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Mean elevation")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
elev_well_plot ## no change in mean elevation between TP1 and TP2
ggsave(elev_well_plot, file="Graphs/Mean_elevation_shift_well_hecs.png")

length(which(mean_elev_well2$elev_shift>0)) ## 39 species moving uphill
length(which(mean_elev_well2$elev_shift<0)) ## 26 species moving downhill
## categorise species as either moving uphill or downhill
mean_elev_well2$shift_direction <- ifelse(mean_elev_well2$elev_shift>0, "Uphill", "Downhill")

## add in elevation results to direction and distance shifts
centroid_elev_well <- merge(mean_elev_well2, range_centroid_well2, by="Common_Name")

## plot result coloured by elevational shift
well_centroid2 <- ggplot(data=centroid_elev_well,
                        aes(x=direction, y=distance, colour=shift_direction)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  ggtitle("Well-recorded hectads") +
  theme_bw()
well_centroid2
ggsave(well_centroid2, file="Graphs/multidirectional_elev_shifts_well.png")


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
tp1_species_percentiles <- nmrsdata_well %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
  summarise(temp_2.5_perc = quantile(temperature, probs = 0.025),
            temp_97.5_perc = quantile(temperature, probs = 0.975))
nmrsdata_well <- merge(nmrsdata_well, tp1_species_percentiles, by="Common_Name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_well %>% group_by(Common_Name, time_period) %>%
  filter(temperature >= temp_2.5_perc & temperature <= temp_97.5_perc)

## calculate centroids of each
temp_centroid_well <- filtered %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
temp_centroid_well2 <- temp_centroid_well %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2

## remove NAs
temp_centroid_well2 <- na.omit(temp_centroid_well2) ## 54 species

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
temp_centroid_well2$distance <- sqrt((temp_centroid_well2$mean_northing_TP1-temp_centroid_well2$mean_northing_TP2)^2 + 
                                      (temp_centroid_well2$mean_easting_TP1-temp_centroid_well2$mean_easting_TP2)^2) ## pythagoras equation
temp_centroid_well2$distance <- temp_centroid_well2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
temp_centroid_well2$bearing <- bearing(temp_centroid_well2[,c(6,4)], temp_centroid_well2[,c(7,5)])
temp_centroid_well2$direction <- (temp_centroid_well2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
temp_centroid_well2$radians <- NISTdegTOradian(temp_centroid_well2$direction)
x <- deg(circ.mean(temp_centroid_well2$radians)) # -36, add full circle
(x+360)%%360 ## 348
# median & quartiles for distances
median(temp_centroid_well2$distance) ## 56.21m 
quantile(temp_centroid_well2$distance) ## 25th = 33.67km, 75th = 95.6km

## plot result
ybreaks <- c(0,50,100,150,200,250,300)
well_temp_centroid <- ggplot(data=temp_centroid_well2,
                            aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) +
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F,
            size = 3) +
  theme_bw() +
  theme( panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
         axis.title.x = element_blank(), axis.title.y = element_blank())
well_temp_centroid
ggsave(well_temp_centroid, file="Graphs/multidirectional_temp_well.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
temp_radians <- temp_centroid_well2[,c("Common_Name", "radians")]
colnames(temp_radians)[2] <- "temp_radians"
range_radians <- range_centroid_well2[,c("Common_Name", "radians")]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(temp_radians, range_radians, by="Common_Name", all.x=TRUE)
circ.cor(radians$temp_radians, radians$range_radians, test=TRUE)
## NON-SIGNIFICANT - species are NOT moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.22

###################
# 2b: Mean summer temperature
###################

#### Direction and distance of shift based on mean annual temperature
## 2.5th and 97.5th percentiles for each species in TP1
detach(package:Rmisc)
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_well %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
  summarise(summer_temp_2.5_perc = quantile(summer_temperature, probs = 0.025),
            summer_temp_97.5_perc = quantile(summer_temperature, probs = 0.975))
nmrsdata_well <- merge(nmrsdata_well, tp1_species_percentiles, by="Common_Name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_well %>% group_by(Common_Name, time_period) %>%
  filter(summer_temperature >= summer_temp_2.5_perc & summer_temperature <= summer_temp_97.5_perc)

## calculate centroids of each
summer_temp_centroid_well <- filtered %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
summer_temp_centroid_well2 <- summer_temp_centroid_well %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2

## remove NAs
summer_temp_centroid_well2 <- na.omit(summer_temp_centroid_well2) ## 56 species

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
summer_temp_centroid_well2$distance <- sqrt((summer_temp_centroid_well2$mean_northing_TP1-summer_temp_centroid_well2$mean_northing_TP2)^2 + 
                                             (summer_temp_centroid_well2$mean_easting_TP1-summer_temp_centroid_well2$mean_easting_TP2)^2) ## pythagoras equation
summer_temp_centroid_well2$distance <- summer_temp_centroid_well2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
summer_temp_centroid_well2$bearing <- bearing(summer_temp_centroid_well2[,c(6,4)], summer_temp_centroid_well2[,c(7,5)])
summer_temp_centroid_well2$direction <- (summer_temp_centroid_well2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
summer_temp_centroid_well2$radians <- NISTdegTOradian(summer_temp_centroid_well2$direction)
x <- deg(circ.mean(summer_temp_centroid_well2$radians)) # -36, add full circle
(x+360)%%360 ## 329
# median & quartiles for distances
median(summer_temp_centroid_well2$distance) ## 43.7km 
quantile(summer_temp_centroid_well2$distance) ## 25th = 25.1km, 75th = 75.8km

## plot result
ybreaks <- c(0,50,100,150,200,250,300)
well_summer_temp_centroid <- ggplot(data=summer_temp_centroid_well2,
                                   aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F,
            size = 3) +
  theme_bw() +
  theme( panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
         axis.title.x = element_blank(), axis.title.y = element_blank())
well_summer_temp_centroid
ggsave(well_summer_temp_centroid, file="Graphs/multidirectional_summer_temp_well.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
temp_radians <- summer_temp_centroid_well2[,c("Common_Name", "radians")]
colnames(temp_radians)[2] <- "temp_radians"
range_radians <- range_centroid_well2[,c("Common_Name", "radians")]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(temp_radians, range_radians, by="Common_Name", all.x=TRUE)
circ.cor(radians$temp_radians, radians$range_radians, test=TRUE)
## significant - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.68
## very similar to mean annual temperature

# 2c: Mean total precipitation
#### Direction and distance of shift based on mean annual temperature
## 2.5th and 97.5th percentiles for each species in TP1
detach(package:Rmisc)
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_well %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
  summarise(precip_2.5_perc = quantile(total_precip, probs = 0.025),
            precip_97.5_perc = quantile(total_precip, probs = 0.975))
nmrsdata_well <- merge(nmrsdata_well, tp1_species_percentiles, by="Common_Name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_well %>% group_by(Common_Name, time_period) %>%
  filter(total_precip >= precip_2.5_perc & total_precip <= precip_97.5_perc)

## calculate centroids of each
precip_centroid_well <- filtered %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
precip_centroid_well2 <- precip_centroid_well %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2

## remove NAs
precip_centroid_well2 <- na.omit(precip_centroid_well2) ## 67 species (5 removed)

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
precip_centroid_well2$distance <- sqrt((precip_centroid_well2$mean_northing_TP1-precip_centroid_well2$mean_northing_TP2)^2 + 
                                        (precip_centroid_well2$mean_easting_TP1-precip_centroid_well2$mean_easting_TP2)^2) ## pythagoras equation
precip_centroid_well2$distance <- precip_centroid_well2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
precip_centroid_well2$bearing <- bearing(precip_centroid_well2[,c(6,4)], precip_centroid_well2[,c(7,5)])
precip_centroid_well2$direction <- (precip_centroid_well2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
precip_centroid_well2$radians <- NISTdegTOradian(precip_centroid_well2$direction)
x <- deg(circ.mean(precip_centroid_well2$radians)) # -36, add full circle
(x+360)%%360 ## 334
# median & quartiles for distances
median(precip_centroid_well2$distance) ## 39.06km 
quantile(precip_centroid_well2$distance) ## 25th = 23.7km, 75th = 76.01km

## plot result
ybreaks <- c(0,50,100,150,200,250,300)
well_precip_centroid <- ggplot(data=precip_centroid_well2,
                              aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  coord_polar() +
  scale_y_continuous(breaks = ybreaks) + 
  geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F,
            size = 3) +
  theme_bw() +
  theme( panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
         axis.title.x = element_blank(), axis.title.y = element_blank())
well_precip_centroid
ggsave(well_precip_centroid, file="Graphs/multidirectional_total_precip_well.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
precip_radians <- precip_centroid_well2[,c("Common_Name", "radians")]
colnames(precip_radians)[2] <- "precip_radians"
range_radians <- range_centroid_well2[,c("Common_Name", "radians")]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(precip_radians, range_radians, by="Common_Name", all.x=TRUE)
circ.cor(radians$precip_radians, radians$range_radians, test=TRUE)
## significant - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.73

## put all circular plots together
library(ggpubr)
circular_well <- ggarrange(well_centroid, well_temp_centroid, well_summer_temp_centroid, well_precip_centroid,
                          labels = c("(a)", "(b)", "(c)", "(d)"),
                          ncol = 2, nrow = 2)
ggsave(circular_well, file="Graphs/circlar_range_temp_summer_precip_well.png", height=12, width=10)

#### 3. Define colonised, persisted, and extirpated hectads
nmrsdata_well$Recorded <- 1
nmrsdata_well_expand <- nmrsdata_well %>% expand(Common_Name, time_period, Hectad)
recorded <- nmrsdata_well[,c("Common_Name", "time_period", "Hectad", "Recorded")]
nmrsdata_well_expand <- merge(nmrsdata_well_expand, recorded, by=c("Common_Name", "time_period", "Hectad"), all.x=TRUE)
nmrsdata_well_expand[is.na(nmrsdata_well_expand)] <- 0 ## change NAs to zero == species was NOT recorded at this hectad in this time period
## change to long format
nmrsdata_well_expand2 <- nmrsdata_well_expand %>%
  spread(time_period, Recorded) ## each species has 1424 rows = the number of recorded hectads
## first remove rows where TP1 AND TP2 == 0 (this a hectad where a species was never recorded)
nmrsdata_well_expand2<-nmrsdata_well_expand2[!(nmrsdata_well_expand2$TP1==0 & nmrsdata_well_expand2$TP2==0),]
nmrsdata_well_expand2$Hectad_category <- case_when(
  nmrsdata_well_expand2$TP1==0 & nmrsdata_well_expand2$TP2==1 ~ "Colonisation",
  nmrsdata_well_expand2$TP1==1 & nmrsdata_well_expand2$TP2==0 ~ "Extirpation",
  TRUE ~ "Persistence"
)

## plot results on a map
## add in lat/lon info
lat_lon <- unique(nmrsdata_well[,c("Hectad", "lat", "lon","elevation10x10km","elevation10x10km_SD")])
nmrsdata_well_expand2 <- merge(nmrsdata_well_expand2, lat_lon, by="Hectad", all.x=TRUE)
worldmap = map_data('world')
ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  theme_void() +
  geom_point(data =nmrsdata_well_expand2,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat), colour=Hectad_category), size=2) +
  theme(title = element_text(size = 12))

## find proportion of species at each extirpated hectad - where are extirpations occurring? 
extirpated_hecs1 <- nmrsdata_well_expand2 %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
  dplyr::summarise(tot_sp = n())
extirpated_hecs2 <- nmrsdata_well_expand2 %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
  filter(Hectad_category=="Extirpation") %>%
  dplyr::summarise(extir_sp = n())
extirpated_hecs <- merge(extirpated_hecs1, extirpated_hecs2, by=c("Hectad", "lat", "lon", "elevation10x10km", "elevation10x10km_SD"))
extirpated_hecs$extir_prop <- extirpated_hecs$extir_sp/extirpated_hecs$tot_sp
## very long way of doing this - gave up on finding a quicker way!

## heatmap
extir_hecs <- ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  theme_void() +
  geom_point(data =extirpated_hecs,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat), colour=extir_prop), size=1) +
  scale_color_viridis_c(name="Proportion of species") + 
  theme(title = element_text(size = 12))
## extirpated hectads with the most species lost tend to be in Northern Scotland/Pennines
extir_hecs
ggsave(extir_hecs, file="Maps/Well_extirpated_hecs.png")


#################################################################
#### How does the climate and elevation of extirpated hectads compared to all recorded hectads?

## proportion plots that Andy suggested
well_hecs <- nmrsdata_well[,c("Hectad","elevation10x10km","elevation10x10km_SD")] ## all rec hecs: 1424
well_hecs <- well_hecs %>% distinct(Hectad, .keep_all = TRUE) ## different way of unique
## subset 'well' extirpated hecs which have lost at least 25% of species
well_extir_hecs <- extirpated_hecs[extirpated_hecs$extir_prop>=0.25,] ## 504 hectads

## add in climate data to each data frame
temp <- read.csv("Data/NMRS/All_NMRS_hectads_annual_temperature.csv", header=TRUE)
summer_temp <- read.csv("Data/NMRS/All_NMRS_hectads_summer_temperature.csv", header=TRUE)
precip <- read.csv("Data/NMRS/All_NMRS_hectads_total_precipitation.csv", header=TRUE)
## put these together
df_list <- list(temp, summer_temp, precip)
climate <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)
## only want TP2 climate data
climate_tp2 <- climate[climate$time_period=="TP2",]

## merge in with all recorded hectads and then all well extirpated hectads
well_hecs <- merge(well_hecs, climate_tp2, by="Hectad", all.x=)
well_extir_hecs <- merge(well_extir_hecs, climate_tp2, by="Hectad", all.x=)

#### Annual temperature
# split temperature into 10 bins
well_hecs$temperature_cat <- cut(well_hecs$temperature, breaks = 10, labels=FALSE)

well_hecs_temp_prop <- well_hecs %>% group_by(temperature_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(well_hecs_temp_prop$freq ~ well_hecs_temp_prop$temperature_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(temperature_cat = cut(well_hecs$temperature, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(temperature_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)
min_max$cat <- 1:10

medians <- well_hecs %>% group_by(temperature_cat) %>% dplyr::summarise(median=median(temperature))

extirpated_hecs2 <- well_extir_hecs %>% mutate(temperature_cat = case_when(
  between(temperature, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(temperature, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(temperature, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(temperature, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(temperature, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(temperature, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(temperature, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(temperature, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(temperature, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(temperature_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$temperature_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("temperature_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
well_hecs_temp_prop <- well_hecs_temp_prop[,c("temperature_cat", "freq")]
well_hecs_temp_prop$hec_cat <- "All hectads"

temp_prop <- rbind(extirpated_hecs_prop, well_hecs_temp_prop)
temp_prop <- merge(temp_prop, medians, by="temperature_cat")
## plot line graph
well_extir_temp <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Mean annual temperature", y="Proportion of hectads") +
  scale_x_continuous(limits=c(4,12)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
well_extir_temp
ggsave(well_extir_temp, file="Graphs/Well_hecs_extirpated_temperature.png")


#### Summer temperature
# split temperature into 10 bins
well_hecs$summer_temperature_cat <- cut(well_hecs$summer_temperature, breaks = 10, labels=FALSE)

well_hecs_temp_prop <- well_hecs %>% group_by(summer_temperature_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(well_hecs_temp_prop$freq ~ well_hecs_temp_prop$summer_temperature_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(summer_temperature_cat = cut(well_hecs$summer_temperature, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(summer_temperature_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)

medians <- well_hecs %>% group_by(summer_temperature_cat) %>% dplyr::summarise(median=median(summer_temperature))

extirpated_hecs2 <- well_extir_hecs %>% mutate(summer_temperature_cat = case_when(
  between(summer_temperature, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(summer_temperature, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(summer_temperature, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(summer_temperature, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(summer_temperature, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(summer_temperature, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(summer_temperature, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(summer_temperature, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(summer_temperature, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(summer_temperature_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$summer_temperature_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("summer_temperature_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
well_hecs_temp_prop <- well_hecs_temp_prop[,c("summer_temperature_cat", "freq")]
well_hecs_temp_prop$hec_cat <- "All hectads"

temp_prop <- rbind(extirpated_hecs_prop, well_hecs_temp_prop)
temp_prop <- merge(temp_prop, medians, by="summer_temperature_cat")
## plot line graph
well_extir_summer_temp <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Mean annual summer temperature", y="Proportion of hectads") +
  scale_x_continuous(limits=c(10,18)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
well_extir_summer_temp
ggsave(well_extir_summer_temp, file="Graphs/Well_hecs_extirpated_summer_temperature.png")



#### Total precipitation
# split precipitatino into 10 bins
well_hecs$precip_cat <- cut(well_hecs$total_precip, breaks = 10, labels=FALSE)

well_hecs_precip_prop <- well_hecs %>% group_by(precip_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(well_hecs_precip_prop$freq ~ well_hecs_precip_prop$precip_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(precip_cat = cut(well_hecs$total_precip, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(precip_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)

medians <- well_hecs %>% group_by(precip_cat) %>% dplyr::summarise(median=median(total_precip))

extirpated_hecs2 <- well_extir_hecs %>% mutate(precip_cat = case_when(
  between(total_precip, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(total_precip, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(total_precip, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(total_precip, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(total_precip, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(total_precip, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(total_precip, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(total_precip, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(total_precip, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(precip_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$precip_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("precip_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
well_hecs_precip_prop <- well_hecs_precip_prop[,c("precip_cat", "freq")]
well_hecs_precip_prop$hec_cat <- "All hectads"

precip_prop <- rbind(extirpated_hecs_prop, well_hecs_precip_prop)
precip_prop <- merge(precip_prop, medians, by="precip_cat")
## plot line graph
well_extir_precip <- ggplot(precip_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Mean total precipitation (mm)", y="Proportion of hectads") +
  scale_x_continuous(breaks=seq(500,3500, by=500)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
well_extir_precip
ggsave(well_extir_precip, file="Graphs/Well_hecs_extirpated_total_precipitation.png")


#### Elevation
# split elevation into 10 bins
well_hecs$elev_cat <- cut(well_hecs$elevation10x10km, breaks = 10, labels=FALSE)

well_hecs_elev_prop <- well_hecs %>% group_by(elev_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(well_hecs_elev_prop$freq ~ well_hecs_elev_prop$elev_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(elev_cat = cut(well_hecs$elevation10x10km, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(elev_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)

medians <- well_hecs %>% group_by(elev_cat) %>% dplyr::summarise(median=median(elevation10x10km))

extirpated_hecs2 <- well_extir_hecs %>% mutate(elev_cat = case_when(
  between(elevation10x10km, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(elevation10x10km, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(elevation10x10km, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(elevation10x10km, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(elevation10x10km, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(elevation10x10km, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(elevation10x10km, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(elevation10x10km, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(elevation10x10km, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(elev_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$elev_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("elev_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
well_hecs_elev_prop <- well_hecs_elev_prop[,c("elev_cat", "freq")]
well_hecs_elev_prop$hec_cat <- "All hectads"

elev_prop <- rbind(extirpated_hecs_prop, well_hecs_elev_prop)
elev_prop <- merge(elev_prop, medians, by="elev_cat")
## plot line graph
well_extir_elev <- ggplot(elev_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Mean elevation (m)", y="Proportion of hectads") +
  scale_x_continuous(breaks=seq(0,2000, by=200)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
well_extir_elev
ggsave(well_extir_elev, file="Graphs/Well_hecs_extirpated_elevation.png")


#### SD Elevation
# split elevation into 10 bins
well_hecs$elev_cat <- cut(well_hecs$elevation10x10km_SD, breaks = 10, labels=FALSE)

well_hecs_elev_prop <- well_hecs %>% group_by(elev_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(well_hecs_elev_prop$freq ~ well_hecs_elev_prop$elev_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(elev_cat = cut(well_hecs$elevation10x10km_SD, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(elev_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)

medians <- well_hecs %>% group_by(elev_cat) %>% dplyr::summarise(median=median(elevation10x10km_SD))

extirpated_hecs2 <- well_extir_hecs %>% mutate(elev_cat = case_when(
  between(elevation10x10km_SD, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(elevation10x10km_SD, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(elevation10x10km_SD, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(elevation10x10km_SD, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(elevation10x10km_SD, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(elevation10x10km_SD, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(elevation10x10km_SD, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(elevation10x10km_SD, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(elevation10x10km_SD, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(elev_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$elev_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("elev_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
well_hecs_elev_prop <- well_hecs_elev_prop[,c("elev_cat", "freq")]
well_hecs_elev_prop$hec_cat <- "All hectads"

elev_prop <- rbind(extirpated_hecs_prop, well_hecs_elev_prop)
elev_prop <- merge(elev_prop, medians, by="elev_cat")
## plot line graph
well_extir_elev_sd <- ggplot(elev_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Standard deviation elevation (m)", y="Proportion of hectads") +
  scale_x_continuous(limits=c(0,400)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
well_extir_elev_sd
ggsave(well_extir_elev_sd, file="Graphs/Well_hecs_extirpated_elevation_SD.png")


### Are hectads which have warmed more over time more likely to be extirpated?
## calculate difference in climate (tp2 - tp1)
temp_df <- climate[,c("Hectad", "time_period", "temperature")]
climate_diff <- temp_df %>%
  spread(time_period, temperature) %>% 
  mutate(temp_diff = TP2-TP1)

## merge in with all recorded hectads and then all well extirpated hectads
well_hecs <- merge(well_hecs, climate_diff, by="Hectad", all.x=TRUE)
well_extir_hecs <- merge(well_extir_hecs, climate_diff, by="Hectad", all.x=TRUE)

#### Annual temperature
# split temperature into 10 bins
well_hecs$temp_diff_cat <- cut(well_hecs$temp_diff, breaks = 10, labels=FALSE)

well_hecs_temp_prop <- well_hecs %>% group_by(temp_diff_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(well_hecs_temp_prop$freq ~ well_hecs_temp_prop$temp_diff_cat) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads

min_max <- tibble(temp_diff_cat = cut(well_hecs$temp_diff, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(temp_diff_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)
min_max <- arrange(min_max, min)

medians <- well_hecs %>% group_by(temp_diff_cat) %>% dplyr::summarise(median=median(temp_diff))

extirpated_hecs2 <- well_extir_hecs %>% mutate(temp_diff_cat = case_when(
  between(temp_diff, min_max[1, "min"], min_max[1, "max"]) ~ 1,
  between(temp_diff, min_max[2, "min"], min_max[2, "max"]) ~ 2,
  between(temp_diff, min_max[3, "min"], min_max[3, "max"]) ~ 3,
  between(temp_diff, min_max[4, "min"], min_max[4, "max"]) ~ 4,
  between(temp_diff, min_max[5, "min"], min_max[5, "max"]) ~ 5,
  between(temp_diff, min_max[6, "min"], min_max[6, "max"]) ~ 6,
  between(temp_diff, min_max[7, "min"], min_max[7, "max"]) ~ 7,
  between(temp_diff, min_max[8, "min"], min_max[8, "max"]) ~ 8,
  between(temp_diff, min_max[9, "min"], min_max[9, "max"]) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(temp_diff_cat) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$temp_diff_cat) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c("temp_diff_cat", "freq")]
extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
well_hecs_temp_prop <- well_hecs_temp_prop[,c("temp_diff_cat", "freq")]
well_hecs_temp_prop$hec_cat <- "All hectads"

temp_prop <- rbind(extirpated_hecs_prop, well_hecs_temp_prop)
temp_prop <- merge(temp_prop, medians, by="temp_diff_cat")
## plot line graph
well_extir_temp_diff <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat), lwd=1) + 
  labs(x="Difference in temperature", y="Proportion of hectads") +
  scale_x_continuous(limits=c(0.25,1.15)) +
  scale_y_continuous(breaks=seq(0,1, by=0.05)) +
  theme_classic() +
  theme(legend.title = element_blank())
well_extir_temp_diff
ggsave(well_extir_temp_diff, file="Graphs/Well_hecs_extirpated_temperature_difference.png")


































# #############################################
# ## direction of centroid shift using minimum convex polygon
# 
# nmrsdata_rec_hecs_tp1 <- nmrsdata_rec[nmrsdata_rec$time_period=="TP1",]
# nmrsdata_rec_hecs_tp2 <- nmrsdata_rec[nmrsdata_rec$time_period=="TP2",]
# 
# ## first time period
# species <- nmrsdata_rec_hecs_tp1$Common_Name
# long <- nmrsdata_rec_hecs_tp1$lon
# lat <- nmrsdata_rec_hecs_tp1$lat
# data <- as.data.frame(cbind(long,lat,species))
# 
# data$lat <- as.numeric(data$lat)
# data$long <- as.numeric(data$long)
# 
# library(letsR)
# PAM <- lets.presab.points(data[,1:2], data[,3], xmn=min(data$long), xmx=max(data$long),
#                           ymn=min(data$lat), ymx=max(data$lat))
# plot(PAM)
# plot(PAM, name="Light Knot Grass")
# summary(PAM)
# centroids_tp1_rec <- lets.midpoint(PAM, planar=FALSE, method="MCC")
# colnames(centroids_tp1_rec) <- c("species","long","lat")
# centroids_tp1_rec$cent <- "yes"
# data$cent <- "no"
# 
# final_centroid <- rbind(data,centroids_tp1_rec)
# species <- unique(final_centroid$species)
# ## plot map for each species with red dot signifying centroid
# ## colour=cent 
# ## save graphs to look at
# worldmap = map_data('world')
# 
# for(i in species) {
#   print(i)
#   # Printing ggplot within for-loop
#   
#   temp_plot <- ggplot() + 
#     geom_polygon(data = worldmap, 
#                  aes(x = long, y = lat, group = group), 
#                  fill = 'gray90', color = 'black') + 
#     coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#     theme_void() + 
#     geom_point(data = final_centroid[final_centroid$species==i,], 
#                aes(x = as.numeric(long), 
#                    y = as.numeric(lat), colour=cent), size=2) +
#     scale_color_manual(values=c("black", "red")) +
#     labs(colour = "Centroid value") 
#   ggtitle(i) +
#     theme(title = element_text(size = 12))
#   
#   ggsave(temp_plot, file=paste0("../../Maps/Centroids/NMRS_centroid_TP1_", i,".png"), width = 20, height = 25, units = "cm")
#   Sys.sleep(2)
# }
# ## all look reasonable - go with this for now
# ## once we've got direction for each species, plot centroid for each TP and check direction
# ## use geom_line to connect the two dots
# 
# ## second time period
# species <- nmrsdata_rec_hecs_tp2$Common_Name
# long <- nmrsdata_rec_hecs_tp2$lon
# lat <- nmrsdata_rec_hecs_tp2$lat
# data <- as.data.frame(cbind(long,lat,species))
# 
# data$lat <- as.numeric(data$lat)
# data$long <- as.numeric(data$long)
# 
# library(letsR)
# PAM <- lets.presab.points(data[,1:2], data[,3], xmn=min(data$long), xmx=max(data$long),
#                           ymn=min(data$lat), ymx=max(data$lat))
# plot(PAM)
# plot(PAM, name="Light Knot Grass")
# summary(PAM)
# centroids_tp2_rec <- lets.midpoint(PAM, planar=FALSE, method="MCC")
# colnames(centroids_tp2_rec) <- c("species","long","lat")
# centroids_tp2_rec$cent <- "yes"
# data$cent <- "no"
# 
# final_centroid <- rbind(data,centroids_tp2_rec)
# species <- unique(final_centroid$species)
# ## plot map for each species with red dot signifying centroid
# ## colour=cent 
# ## save graphs to look at
# worldmap = map_data('world')
# 
# for(i in species) {
#   print(i)
#   # Printing ggplot within for-loop
#   
#   temp_plot <- ggplot() + 
#     geom_polygon(data = worldmap, 
#                  aes(x = long, y = lat, group = group), 
#                  fill = 'gray90', color = 'black') + 
#     coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#     theme_void() + 
#     geom_point(data = final_centroid[final_centroid$species==i,], 
#                aes(x = as.numeric(long), 
#                    y = as.numeric(lat), colour=cent), size=2) +
#     scale_color_manual(values=c("black", "red")) +
#     labs(colour = "Centroid value") 
#   ggtitle(i) +
#     theme(title = element_text(size = 12))
#   
#   ggsave(temp_plot, file=paste0("../../Maps/Centroids/NMRS_centroid_TP2_", i,".png"), width = 20, height = 25, units = "cm")
#   Sys.sleep(2)
# }
# 
# ## remove cent and species column
# species <- centroids_tp1_rec$species
# centroids_tp1_rec <- centroids_tp1_rec[,2:3]
# centroids_tp2_rec <- centroids_tp2_rec[,2:3]
# 
# 
# ## measure distance? and direction of shift
# library(geosphere)
# x <- bearing(centroids_tp1_rec, centroids_tp2_rec)
# ## add in time period
# centroids_tp1_rec$time_period <- "TP1"
# centroids_tp2_rec$time_period <- "TP2"
# 
# centroids_rec_hecs <- rbind(centroids_tp1_rec, centroids_tp2_rec)
# centroids_rec_hecs$species <- species
# 
# direction <- data.frame(bearing=x, species=species)
# course <- (x + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
# direction$course  <- course
# ## reorder columns
# direction <- direction[,c(2,1,3)]
# 
# ## plot each species centroid in TP1 and TP2 to check direction
# ## save graphs to look at
# worldmap = map_data('world')
# 
# for(i in species) {
#   print(i)
#   # Printing ggplot within for-loop
#   
#   temp_plot <- ggplot() + 
#     geom_polygon(data = worldmap, 
#                  aes(x = long, y = lat, group = group), 
#                  fill = 'gray90', color = 'black') + 
#     coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#     theme_void() + 
#     geom_point(data = centroids_rec_hecs[centroids_rec_hecs$species==i,], 
#                aes(x = as.numeric(long), 
#                    y = as.numeric(lat), colour=time_period), size=2) +
#     geom_line() +
#     labs(colour = "Time period") 
#   ggtitle(i) +
#     theme(title = element_text(size = 12))
#   
#   ggsave(temp_plot, file=paste0("../../Maps/Centroids/NMRS_centroid_shift_", i,".png"), width = 20, height = 25, units = "cm")
#   Sys.sleep(2)
# }
# 
# ## get each species classed into 8-point compass directions
# ## N = 337.5 - 22.5
# ## NE = 22.5 - 67.5
# ## E = 67.5 - 112.5
# ## SE = 112.5 - 157.5
# ## S = 157.5 - 202.5
# ## SW = 202.5 - 247.5
# ## W = 247.5 - 292.5
# ## NW = 292.5 - 337.5
# 
# ## long way of doing it but it works
# direction$compass_direction <- NULL
# NE <- direction %>% filter(between(course,22.5,67.5) )
# NE$compass_direction <- "NE"
# E <- direction %>% filter(between(course,67.5,112.5) )
# E$compass_direction <- "E"
# SE <- direction %>% filter(between(course,112.5,157.5) )
# SE$compass_direction <- "SE"
# S <- direction %>% filter(between(course,157.5,202.5) )
# S$compass_direction <- "S"
# SW <- direction %>% filter(between(course,202.5,247.5) )
# SW$compass_direction <- "SW"
# W <- direction %>% filter(between(course,247.5,292.5) )
# W$compass_direction <- "W"
# NW <- direction %>% filter(between(course,292.5,337.5) )
# NW$compass_direction <- "NW"
# N <- direction %>% filter(course>337.5 | course<22.5) 
# N$compass_direction <- "N"
# 
# directions_rec_hecs <- rbind(NE,E,SE,S,SW,W,NW,N)
# ## nothing has moved in east direction
# #############################################

