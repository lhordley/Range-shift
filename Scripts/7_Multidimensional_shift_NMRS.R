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
rec_centroid <- ggplot(data=range_centroid_rec2,
                       aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  theme_bw()
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
rec_temp_centroid <- ggplot(data=temp_centroid_rec2,
                       aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  theme_bw()
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
temp_centroid_rec2 <- na.omit(temp_centroid_rec2) ## 66 species (7 removed)

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
(x+360)%%360 ## 327
# median & quartiles for distances
median(temp_centroid_rec2$distance) ## 46.3km 
quantile(temp_centroid_rec2$distance) ## 25th = 27.18km, 75th = 72.8km

## plot result
rec_summer_temp_centroid <- ggplot(data=temp_centroid_rec2,
                            aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  theme_bw()
rec_summer_temp_centroid
ggsave(rec_summer_temp_centroid, file="Graphs/multidirectional_summer_temp_rec.png")

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
# median & quartiles for distances
median(precip_centroid_rec2$distance) ## 35.26km 
quantile(precip_centroid_rec2$distance) ## 25th = 20.12km, 75th = 57.12km

## plot result
rec_precip_centroid <- ggplot(data=precip_centroid_rec2,
                                   aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  theme_bw()
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


#### 3. Define colonised, persisted, and extirpated hectads
nmrsdata_rec$Recorded <- 1
nmrsdata_rec_expand <- nmrsdata_rec %>% expand(Common_Name, time_period, Hectad)
recorded <- nmrsdata_rec[,c(1:3,25)]
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
lat_lon <- unique(nmrsdata_rec[,c(2,6:9)])
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
extirpated_hecs <- nmrsdata_rec_expand2 %>% group_by(Hectad) %>%
                    dplyr::summarise(cnt = n()) %>%
                    filter(Hectad_category=="Extirpation")
                    dplyr::summarise(freq = n() / cnt)
                    ###################################
                    ######################################
                    ######################################
                    ######################################
                    ########################################
                    #######################################
                    ###################################
extirpated_hecs <- nmrsdata_rec_expand2[nmrsdata_rec_expand2$Hectad_category=="Extirpation",]
extirpated_hecs2 <- extirpated_hecs %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
  dplyr::summarise(n_sp=n()) ## 841 hectads
## heatmap
extir_hecs <- ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  theme_void() +
  geom_point(data =extirpated_hecs2,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat), colour=n_sp), size=1) +
  scale_color_viridis_c(name="Number of species") + 
  theme(title = element_text(size = 12))
## extirpated hectads with the most species lost tend to be in Northern Scotland/Pennines
extir_hecs
ggsave(extir_hecs, file="Maps/Rec_extirpated_hecs.png")

## proportion plots that Andy suggested
rec_hecs <- unique(nmrsdata_rec[,c(2,8)]) ## all well hecs: 1424
rec_hecs$elev_cat <- cut(rec_hecs$elevation10x10km, breaks = 10)
rec_hecs$elev_cat_int <- cut(rec_hecs$elevation10x10km, breaks = 10, labels=FALSE)

rec_hecs_prop <- rec_hecs %>% group_by(elev_cat_int) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(rec_hecs_prop$freq ~ rec_hecs_prop$elev_cat_int) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads
## only take those with the highers # of extirpations (10 or more species per hectad?)
extirpated_hecs2 <- extirpated_hecs2[extirpated_hecs2$n_sp>=10,] ## 93 hectads

min_max <- tibble(elev_cat = cut(rec_hecs$elevation10x10km, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(elev_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)

extirpated_hecs2 <- extirpated_hecs2 %>% mutate(elev_cat_int = case_when(
  between(elevation10x10km, 0, 138) ~ 1,
  between(elevation10x10km, 138, 276) ~ 2,
  between(elevation10x10km, 276, 414) ~ 3,
  between(elevation10x10km, 414, 551) ~ 4,
  between(elevation10x10km, 551, 689) ~ 5,
  between(elevation10x10km, 689, 827) ~ 6,
  between(elevation10x10km, 827, 965) ~ 7,
  between(elevation10x10km, 965, 1100) ~ 8,
  between(elevation10x10km, 1100, 1240) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(elev_cat_int) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$elev_cat_int) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c(1,3)]
extirpated_hecs_prop$hec_cat <- "extirpated"
rec_hecs_prop <- rec_hecs_prop[,c(1,3)]
rec_hecs_prop$hec_cat <- "all"

elevation_prop <- rbind(extirpated_hecs_prop, rec_hecs_prop)
## plot line graph
rec_extir_elev <- ggplot(elevation_prop, aes(x=elev_cat_int, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat)) + 
  labs(x="Elevation categories", y="Proportion of hectads") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  theme_classic()
rec_extir_elev
ggsave(rec_extir_elev, file="Graphs/Rec_hecs_extirpated_elevation.png")

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


###############################################################
## well & heavily recorded hectads

## filter to recorded hectads only (i.e. those recorded once in both time periods)
well_heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded" | RECORDING.LEVEL=="Well recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
well_heavy_hecs <- well_heavy_hecs[which(well_heavy_hecs$n_row > 1), ]
well_heavy_hecs$n_row <- NULL ## 747 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_wh <- nmrsdata[which(nmrsdata$Hectad %in% well_heavy_hecs$HECTAD), ]
length(unique(nmrsdata_wh$Hectad)) ## 639 hectads
length(unique(nmrsdata_wh$Common_Name)) ## 69 species

##### Mean elevation shift
mean_elev_well <- nmrsdata_wh %>% group_by(Common_Name, time_period) %>%
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

mean(mean_elev_well2$elev_shift) ## 42.03 uphill
sd(mean_elev_well2$elev_shift) ## 149.4m 

# Shapiro-Wilk normality test for the differences
shapiro.test(mean_elev_well2$elev_shift) # => p-value = 0.02 - distribution is NOT normal
qqnorm(mean_elev_well2$elev_shift)
qqline(mean_elev_well2$elev_shift) ## looks ok

## one sample t-test to test whether elevation shifts differ significantly from 0
t.test(mean_elev_well2$elev_shift) ## significant
wilcox.test(mean_elev_well2$elev_shift) ## significant

## check for normality 
# compute the difference
d <- with(mean_elev_well, 
          mean_elev[time_period == "TP1"] - mean_elev[time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = <0.001 - distribution is NOT normal
qqnorm(d)
qqline(d) ## not good

## wilcoxon signed-rank test instead
elev_well <- wilcox.test(mean_elev ~ time_period, data = mean_elev_well, paired = TRUE)
elev_well ## significant

## plot result
elev_well_plot <- ggpaired(mean_elev_well, x = "time_period", y = "mean_elev",
                        color = "time_period", line.color = "gray", line.size = 0.4,
                        palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Mean elevation")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
elev_well_plot
ggsave(elev_well_plot, file="Graphs/Mean_elevation_shift_well_hecs.png")


length(which(mean_elev_well2$elev_shift>0)) ## 44 species moving uphill (66%)
length(which(mean_elev_well2$elev_shift<0)) ## 21 species moving downhill (33%)

mean_elev_well2$shift_direction <- ifelse(mean_elev_well2$elev_shift>0, "Uphill", "Downhill")

####### Multidirectional centroid range shift

# take mean easting and northing and lat and lon across each species occupied hectads for each time period = range centroid
# following Gillings et al 2015
range_centroid_wh <- nmrsdata_wh %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))
## change orientation of data
range_centroid_wh2 <- range_centroid_wh %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
## remove NAs - some species are only found in one TP 
range_centroid_wh2 <- na.omit(range_centroid_wh2) ## 67 species
# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
range_centroid_wh2$distance <- sqrt((range_centroid_wh2$mean_northing_TP1-range_centroid_wh2$mean_northing_TP2)^2 + 
                                      (range_centroid_wh2$mean_easting_TP1-range_centroid_wh2$mean_easting_TP2)^2) ## pythagoras equation
range_centroid_wh2$distance <- range_centroid_wh2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
range_centroid_wh2$bearing <- bearing(range_centroid_wh2[,c(6,4)], range_centroid_wh2[,c(7,5)])
range_centroid_wh2$direction <- (range_centroid_wh2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360

## add in elevation results
centroid_elev_well <- merge(mean_elev_well2, range_centroid_wh2, by="Common_Name")

## plot result
wh_centroid <- ggplot(data=centroid_elev_well,
                      aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  theme_bw()
wh_centroid
ggsave(wh_centroid, file="Graphs/multidirectional_shifts_well.png")

## plot result with elevation
wh_centroid2 <- ggplot(data=centroid_elev_well,
                      aes(x=direction, y=distance, colour=shift_direction)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  theme_bw()
wh_centroid2
ggsave(wh_centroid2, file="Graphs/multidirectional_elev_shifts_well.png")

## plot barchart of elevational shifts
centroid_elev_well$effort_level <- "Well recorded hectads"
centroid_elev_rec$effort_level <- "Recorded hectads"

centroid_elev <- rbind(centroid_elev_rec, centroid_elev_well)
elev_sum <- centroid_elev %>% group_by(effort_level, shift_direction) %>% dplyr::summarise(n_sp=n())
elev_sum$tot_sp <- c(72,72,67,67)
elev_sum$perc_change <- (elev_sum$n_sp/elev_sum$tot_sp)*100  

elev_shift <- ggplot(elev_sum, aes(fill=effort_level, y=perc_change, x=shift_direction)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(y="Percentage change", x="Direction of elevational shift") +
  theme_classic()
elev_shift

# use circ.mean from the circstats package to find a mean direction & confidence limits across all species
library(NISTunits)
library(CircStats)
range_centroid_wh2$radians <- NISTdegTOradian(range_centroid_wh2$direction)
x <- deg(circ.mean(range_centroid_wh2$radians)) # -53, add full circle
(x+360)%%360 ## 320
# median & quartiles for distances
median(range_centroid_wh2$distance) ## 37.3km 
quantile(range_centroid_wh2$distance) ## 25th = 22.8km, 75th = 68km

## Rayleigh's test
## if significant, the distribution of directions differ significantly from an even distribution
## if R = 0 (the data are completely spread around the circle)
## if R = 1 (the data are completely concentrated on one point)
r.test(range_centroid_wh2$radians)
## p < 0.001
## r = 0.36

library(circular)

###############
## Define colonised, persisted, and extirpated hectads
## Where have range shifts caused extirpations? 
nmrsdata_wh$Recorded <- 1
nmrsdata_well_expand <- nmrsdata_wh %>% expand(Common_Name, time_period, Hectad)
recorded <- nmrsdata_wh[,c(1,2,11,12)]
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
lat_lon <- unique(nmrsdata_wh[,c(2,6:9)])
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

## find number of species at each extirpated hectad - where are extirpations occurring? 
# hectad, lat, lon, hectad_category, n_sp
extirpated_hecs <- nmrsdata_well_expand2[nmrsdata_well_expand2$Hectad_category=="Extirpation",]
extirpated_hecs2 <- extirpated_hecs %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
  dplyr::summarise(n_sp=n()) ## 496 hectads

## heatmap
extir_hecs <- ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  theme_void() +
  geom_point(data =extirpated_hecs2,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat), colour=n_sp), size=1) +
  scale_color_viridis_c(name="Number of species") + 
  theme(title = element_text(size = 12))
## extirpated hectads with the most species lost tend to be in Northern Scotland
extir_hecs
ggsave(extir_hecs, file="Maps/Well_rec_extirpated_hecs.png")

hist(extirpated_hecs2$elevation10x10km) ## more hectads extirpated at lower elevatin
hist(extirpated_hecs2$elevation10x10km_SD) ## more hectads extirpated with lower variation
## species are moving away from hectads that have fewer choices of where to go
## is this different for persisted hectads?
persist_hecs <- nmrsdata_well_expand2[nmrsdata_well_expand2$Hectad_category=="Persistence",]
persist_hecs2 <- persist_hecs %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
  dplyr::summarise(n_sp=n()) ## 496 hectads
hist(persist_hecs2$elevation10x10km) ## lots of persisted hectads at lower elevation
hist(persist_hecs2$elevation10x10km_SD) ## more hectads persisted with lower variation
## species don't seem to be persisting more in hectads with lower elevation

## now find number of extirpated, colonised and persisted hectads for each species
species_hec_cat <- nmrsdata_well_expand2 %>% group_by(Common_Name, Hectad_category) %>%
  dplyr::summarise(n_hecs=n()) ## 496 hectads
## find number of hectads for each species in TP1 - express the extirpated/colonised/persisted as a % of these
sp_hist_hectads <- nmrsdata_wh %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
  dplyr::summarise(hist_hecs=n())
species_hec_cat <- merge(species_hec_cat, sp_hist_hectads, by="Common_Name", all.x=TRUE)
## change NA to zero (for Northern Dart - not found in TP1)
species_hec_cat <- na.omit(species_hec_cat) ## Remove Northern Dart - only found in one hectad in TP2
species_hec_cat$perc <- (species_hec_cat$n_hecs/species_hec_cat$hist_hecs)*100

boxplot(perc ~ Hectad_category, data=species_hec_cat)

## what elevation are the hecatds that are being extirpated at? Lower hectads more likely to be extirpated or higher?

## proportion plots that Andy suggested
well_hecs <- unique(nmrsdata_wh[,c(2,8)]) ## all well hecs: 639
well_hecs$elev_cat <- cut(well_hecs$elevation10x10km, breaks = 10)
well_hecs$elev_cat_int <- cut(well_hecs$elevation10x10km, breaks = 10, labels=FALSE)

well_hecs_prop <- well_hecs %>% group_by(elev_cat_int) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(well_hecs_prop$freq ~ well_hecs_prop$elev_cat_int) ## lower proportion of hectads with higher elevation
## same again for extirpated hectads
## only take those with the highers # of extirpations (10 or more species per hectad?)
extirpated_hecs2 <- extirpated_hecs2[extirpated_hecs2$n_sp>=10,]

min_max <- tibble(elev_cat = cut(well_hecs$elevation10x10km, breaks = 10)) %>% 
  mutate(x_tmp = str_sub(elev_cat, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double)
min_max <- unique(min_max)

extirpated_hecs2 <- extirpated_hecs2 %>% mutate(elev_cat_int = case_when(
  between(elevation10x10km, 0, 124) ~ 1,
  between(elevation10x10km, 124, 249) ~ 2,
  between(elevation10x10km, 249, 373) ~ 3,
  between(elevation10x10km, 373, 498) ~ 4,
  between(elevation10x10km, 498, 622) ~ 5,
  between(elevation10x10km, 622, 746) ~ 6,
  between(elevation10x10km, 746, 871) ~ 7,
  between(elevation10x10km, 871, 995) ~ 8,
  between(elevation10x10km, 995, 1120) ~ 9,
  TRUE ~ 10
))

extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(elev_cat_int) %>% 
  dplyr::summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt), 3)
plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$elev_cat_int) ## lower proportion of hectads with higher elevation

extirpated_hecs_prop <- extirpated_hecs_prop[,c(1,3)]
extirpated_hecs_prop$hec_cat <- "extirpated"
well_hecs_prop <- well_hecs_prop[,c(1,3)]
well_hecs_prop$hec_cat <- "all"

elevation_prop <- rbind(extirpated_hecs_prop, well_hecs_prop)
## plot line graph
well_extir_elev <- ggplot(elevation_prop, aes(x=elev_cat_int, y=freq, group=hec_cat)) +
  geom_line(aes(linetype=hec_cat)) + 
  labs(x="Elevation categories", y="Proportion of hectads") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  theme_classic()
well_extir_elev
ggsave(well_extir_elev, file="Graphs/Well_hecs_extirpated_elevation.png")

##################### Same for mean annual temperature in TP2
well_hecs <- unique(nmrsdata_wh[,c(2,10:11)]) ## all well hecs: 639
well_hecs <- well_hecs[well_hecs$time_period=="TP2"] ## only 543 hectads
######## NEED TEMPERATURE DATA FOR EXTIRPATED HECTADS


