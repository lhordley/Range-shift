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
hectads <- unique(hec_records$HECTAD) ## 2191
nmrsdata_rec <- nmrsdata[which(nmrsdata$Hectad %in% hectads), ]
length(unique(nmrsdata_rec$Hectad)) ## 1886 hectads
length(unique(nmrsdata_rec$Common_name)) ## 77 species

#### 1. Multidimensional RANGE shift 
# Distance and direction of centroid shift between TP1 and TP2
# take mean easting and northing and lat and lon across each species occupied hectads for each time period = range centroid
# following Gillings et al 2015
range_centroid_rec <- nmrsdata_rec %>% group_by(Common_name, Time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
range_centroid_rec2 <- range_centroid_rec %>%
  gather(key, value, -Common_name, -Time_period) %>%
  unite(col, key, Time_period) %>%
  spread(col, value)

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
range_centroid_rec2$distance <- sqrt((range_centroid_rec2$`mean_northing_1975-1991`-range_centroid_rec2$`mean_northing_2012-2016`)^2 + 
                                       (range_centroid_rec2$`mean_easting_1975-1991`-range_centroid_rec2$`mean_easting_2012-2016`)^2) ## pythagoras equation
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
ggsave(rec_centroid, file="Graphs/multidirectional_shifts_rec_new.png")

# use circ.mean from the circstats package to find a mean direction & confidence limits across all species
library(NISTunits)
library(CircStats)
range_centroid_rec2$radians <- NISTdegTOradian(range_centroid_rec2$direction)
x <- deg(circ.mean(range_centroid_rec2$radians)) # -36, add full circle
(x+360)%%360 ## 325
ci <- vm.bootstrap.ci(range_centroid_rec2$radians) 
deg(c(5.44, 5.89)) ## lower=311.7 upper=337.47

# median & quartiles for distances
median(range_centroid_rec2$distance) ## 40.7km 
quantile(range_centroid_rec2$distance) ## 25th = 23.1km, 75th = 71.9km

## Rayleigh's test
## if significant, the distribution of directions differ significantly from an even distribution
## if R = 0 (the data are completely spread around the circle)
## if R = 1 (the data are completely concentrated on one point)
r.test(range_centroid_rec2$radians)
## p < 0.001
## r = 0.55
## low R = variable directions, but significant means it has significant bias to a certain direction (the mean)

# ## Are species that are moving southwards moving uphill?
# # Test by calculating mean elevation change between time periods
# mean_elev_rec <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>%
#   dplyr::summarise(mean_elev=mean(elevation10x10km))
# ## change orientation
# mean_elev_rec2 <- mean_elev_rec %>%
#   gather(key, value, -Common_Name, -time_period) %>%
#   unite(col, key, time_period) %>%
#   spread(col, value)
# mean_elev_rec2$elev_shift <- mean_elev_rec2$mean_elev_TP2 - mean_elev_rec2$mean_elev_TP1 ## TP2 - TP1
# hist(mean_elev_rec2$elev_shift) ## looks good
# 
# mean(mean_elev_rec2$elev_shift) ## 38.23 m uphill
# sd(mean_elev_rec2$elev_shift) ## 111.27m 
# 
# # Shapiro-Wilk normality test for the differences
# shapiro.test(mean_elev_rec2$elev_shift) # => p-value = 0.02 - distribution is NOT normal
# qqnorm(mean_elev_rec2$elev_shift)
# qqline(mean_elev_rec2$elev_shift) ## not good
# 
# ## Wilcox test to test whether elevation shifts differ significantly from 0
# wilcox.test(mean_elev_rec2$elev_shift) ## significant
# 
# ## wilcoxon signed-rank test - does mean elevation differ between time periods?
# elev_rec <- wilcox.test(mean_elev ~ time_period, data = mean_elev_rec, paired = TRUE)
# elev_rec ## significant
# 
# ## plot result
# elev_rec_plot <- ggpaired(mean_elev_rec, x = "time_period", y = "mean_elev",
#                            color = "time_period", line.color = "gray", line.size = 0.4,
#                            palette = "jco", id="Common_Name")+
#   xlab("Time period")+
#   ylab("Mean elevation")+
#   stat_compare_means(method="wilcox.test", paired = TRUE)
# elev_rec_plot ## higher mean elevation in TP2 compared to TP1
# ggsave(elev_rec_plot, file="Graphs/Mean_elevation_shift_rec_hecs.png")
# 
# length(which(mean_elev_rec2$elev_shift>0)) ## 49 species moving uphill
# length(which(mean_elev_rec2$elev_shift<0)) ## 23 species moving downhill
# ## categorise species as either moving uphill or downhill
# mean_elev_rec2$shift_direction <- ifelse(mean_elev_rec2$elev_shift>0, "Uphill", "Downhill")
# 
# ## add in elevation results to direction and distance shifts
# centroid_elev_rec <- merge(mean_elev_rec2, range_centroid_rec2, by="Common_Name")
# 
# ## plot result coloured by elevational shift
# rec_centroid2 <- ggplot(data=centroid_elev_rec,
#                        aes(x=direction, y=distance, colour=shift_direction)) +
#   geom_segment(aes(xend = direction, yend = 0.1)) +
#   geom_point() +
#   scale_x_continuous(limits = c(0,360),
#                      breaks = seq(0, 360, by = 45),
#                      minor_breaks = seq(0, 360, by = 15)) +
#   # scale_y_log10() +
#   coord_polar() +
#   ggtitle("Recorded hectads") +
#   theme_bw()
# rec_centroid2
# ggsave(rec_centroid2, file="Graphs/multidirectional_elev_shifts_rec.png")


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
tp1_species_percentiles <- nmrsdata_rec %>% group_by(Common_name) %>% filter(Time_period=="1975-1991") %>%
  summarise(temp_2.5_perc = quantile(temperature, probs = 0.025),
            temp_97.5_perc = quantile(temperature, probs = 0.975))
nmrsdata_rec <- merge(nmrsdata_rec, tp1_species_percentiles, by="Common_name")

## find hectads in TP2 which fall within 2.5th and 97.5th percentiles for each species
filtered <- nmrsdata_rec %>% group_by(Common_name, Time_period) %>%
  filter(temperature >= temp_2.5_perc & temperature <= temp_97.5_perc)

## calculate centroids of each
temp_centroid_rec <- filtered %>% group_by(Common_name, Time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
temp_centroid_rec2 <- temp_centroid_rec %>%
  gather(key, value, -Common_name, -Time_period) %>%
  unite(col, key, Time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 2.5th & 97.5th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2
## NAs = slender scotch burnet and silurian
## Grey is also removed - no sites in T2 fall within percentiles

## remove NAs
temp_centroid_rec2 <- na.omit(temp_centroid_rec2) ## 74 species (3 removed)

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
temp_centroid_rec2$distance <- sqrt((temp_centroid_rec2$`mean_northing_1975-1991`-temp_centroid_rec2$`mean_northing_2012-2016`)^2 + 
                                      (temp_centroid_rec2$`mean_easting_1975-1991`-temp_centroid_rec2$`mean_easting_2012-2016`)^2) ## pythagoras equation
temp_centroid_rec2$distance <- temp_centroid_rec2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
temp_centroid_rec2$bearing <- bearing(temp_centroid_rec2[,c(6,4)], temp_centroid_rec2[,c(7,5)])
temp_centroid_rec2$direction <- (temp_centroid_rec2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
temp_centroid_rec2$radians <- NISTdegTOradian(temp_centroid_rec2$direction)
x <- deg(circ.mean(temp_centroid_rec2$radians)) # -36, add full circle
(x+360)%%360 ## 339.95
ci <- vm.bootstrap.ci(temp_centroid_rec2$radians) 
deg(c(-0.54, -0.17)) ## lower=327 upper=347????
# median & quartiles for distances
median(temp_centroid_rec2$distance) ## 69.66km 
quantile(temp_centroid_rec2$distance) ## 25th = 26.79km, 75th = 110.06km

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
ggsave(rec_temp_centroid, file="Graphs/multidirectional_temp_rec_new.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
temp_radians <- temp_centroid_rec2[,c(1,13)]
colnames(temp_radians)[2] <- "temp_radians"
range_radians <- range_centroid_rec2[,c(1,13)]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(temp_radians, range_radians, by="Common_name", all.x=TRUE)
circ.cor(radians$temp_radians, radians$range_radians, test=TRUE)
## significant - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.61

# ###################
# # 2b: Mean summer temperature
# ###################
# 
# #### Direction and distance of shift based on mean annual temperature
# ## 2.5th and 97.5th percentiles for each species in TP1
# detach(package:Rmisc)
# detach(package:plyr)
# tp1_species_percentiles <- nmrsdata_rec %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
#   summarise(summer_temp_2.5_perc = quantile(summer_temperature, probs = 0.025),
#             summer_temp_97.5_perc = quantile(summer_temperature, probs = 0.975))
# nmrsdata_rec <- merge(nmrsdata_rec, tp1_species_percentiles, by="Common_Name")
# 
# ## find hectads in TP2 which fall within 25th and 75th percentiles for each species
# filtered <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>%
#   filter(summer_temperature >= summer_temp_2.5_perc & summer_temperature <= summer_temp_97.5_perc)
# 
# ## calculate centroids of each
# summer_temp_centroid_rec <- filtered %>% group_by(Common_Name, time_period) %>%
#   dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))
# 
# ## change orientation of data
# summer_temp_centroid_rec2 <- summer_temp_centroid_rec %>%
#   gather(key, value, -Common_Name, -time_period) %>%
#   unite(col, key, time_period) %>%
#   spread(col, value)
# ## NAs exist either because:
# # Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# # No hectads fall between 25th and 75th percentiles in TP2
# 
# ## remove NAs
# summer_temp_centroid_rec2 <- na.omit(summer_temp_centroid_rec2) ## 66 species (7 removed)
# 
# # then calculate distance (magnitude) and direction of shift between time periods for each species
# # Use pythagoras equation
# summer_temp_centroid_rec2$distance <- sqrt((summer_temp_centroid_rec2$mean_northing_TP1-summer_temp_centroid_rec2$mean_northing_TP2)^2 + 
#                                       (summer_temp_centroid_rec2$mean_easting_TP1-summer_temp_centroid_rec2$mean_easting_TP2)^2) ## pythagoras equation
# summer_temp_centroid_rec2$distance <- summer_temp_centroid_rec2$distance/1000
# ## this gives distance in km
# 
# # use mean lat and lon values to calculate bearing/direction using geosphere package
# library(geosphere)
# summer_temp_centroid_rec2$bearing <- bearing(summer_temp_centroid_rec2[,c(6,4)], summer_temp_centroid_rec2[,c(7,5)])
# summer_temp_centroid_rec2$direction <- (summer_temp_centroid_rec2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
# summer_temp_centroid_rec2$radians <- NISTdegTOradian(summer_temp_centroid_rec2$direction)
# x <- deg(circ.mean(summer_temp_centroid_rec2$radians)) # -36, add full circle
# (x+360)%%360 ## 328
# ci <- vm.bootstrap.ci(summer_temp_centroid_rec2$radians) 
# deg(c(5.47, 5.93)) ## lower=313 upper=340
# # median & quartiles for distances
# median(summer_temp_centroid_rec2$distance) ## 46.3km 
# quantile(summer_temp_centroid_rec2$distance) ## 25th = 27.18km, 75th = 72.8km
# 
# ## plot result
# ybreaks <- c(0,50,100,150,200,250,300)
# rec_summer_temp_centroid <- ggplot(data=summer_temp_centroid_rec2,
#                             aes(x=direction, y=distance)) +
#   geom_segment(aes(xend = direction, yend = 0.1)) +
#   geom_point() +
#   scale_x_continuous(limits = c(0,360),
#                      breaks = seq(0, 360, by = 45),
#                      minor_breaks = seq(0, 360, by = 15)) +
#   coord_polar() +
#   scale_y_continuous(breaks = ybreaks) + 
#   geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
#             aes(x = x, y = y, label = label),
#             inherit.aes = F,
#             size = 3) +
#   theme_bw() +
#   theme( panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
#          axis.title.x = element_blank(), axis.title.y = element_blank())
# rec_summer_temp_centroid
# ggsave(rec_summer_temp_centroid, file="Graphs/multidirectional_summer_temp_rec.png")
# 
# ## Correlate direction of actual range shift and temperature shift 
# library(CircStats)
# temp_radians <- summer_temp_centroid_rec2[,c(1,13)]
# colnames(temp_radians)[2] <- "temp_radians"
# range_radians <- range_centroid_rec2[,c(1,13)]
# colnames(range_radians)[2] <- "range_radians"
# radians <- merge(temp_radians, range_radians, by="Common_Name", all.x=TRUE)
# circ.cor(radians$temp_radians, radians$range_radians, test=TRUE)
# ## significant - species are moving in a similar direction to what is expected for them
# ## to stay within their temperature envelope
# ## r = 0.752
# ## very similar to mean annual temperature

# 2c: Mean total precipitation
#### Direction and distance of shift based on mean annual temperature
## 2.5th and 97.5th percentiles for each species in TP1
detach(package:Rmisc)
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_rec %>% group_by(Common_name) %>% filter(Time_period=="1975-1991") %>%
  summarise(precip_2.5_perc = quantile(total_precip, probs = 0.025),
            precip_97.5_perc = quantile(total_precip, probs = 0.975))
nmrsdata_rec <- merge(nmrsdata_rec, tp1_species_percentiles, by="Common_name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_rec %>% group_by(Common_name, Time_period) %>%
  filter(total_precip >= precip_2.5_perc & total_precip <= precip_97.5_perc)

## calculate centroids of each
precip_centroid_rec <- filtered %>% group_by(Common_name, Time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))

## change orientation of data
precip_centroid_rec2 <- precip_centroid_rec %>%
  gather(key, value, -Common_name, -Time_period) %>%
  unite(col, key, Time_period) %>%
  spread(col, value)
## NAs exist either because:
# Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# No hectads fall between 25th and 75th percentiles in TP2
# NAs: Grey, Silurian and Scotch Burnet

## remove NAs
precip_centroid_rec2 <- na.omit(precip_centroid_rec2) ## 74 species (3 removed)

# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
precip_centroid_rec2$distance <- sqrt((precip_centroid_rec2$`mean_northing_1975-1991`-precip_centroid_rec2$`mean_northing_2012-2016`)^2 + 
                                        (precip_centroid_rec2$`mean_easting_1975-1991`-precip_centroid_rec2$`mean_easting_2012-2016`)^2) ## pythagoras equation
precip_centroid_rec2$distance <- precip_centroid_rec2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
precip_centroid_rec2$bearing <- bearing(precip_centroid_rec2[,c(6,4)], precip_centroid_rec2[,c(7,5)])
precip_centroid_rec2$direction <- (precip_centroid_rec2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
precip_centroid_rec2$radians <- NISTdegTOradian(precip_centroid_rec2$direction)
x <- deg(circ.mean(precip_centroid_rec2$radians)) # -36, add full circle
(x+360)%%360 ## 329
x.bs <- vm.bootstrap.ci(precip_centroid_rec2$radians)
deg(c(5.52, 5.98)) ## lower=316.27 upper=342.63
library(Directional)
circ.summary(precip_centroid_rec2$radians, rads=TRUE)

r.test(precip_centroid_rec2$radians) ## significant, r=0.578

# median & quartiles for distances
median(precip_centroid_rec2$distance) ## 39.33km 
quantile(precip_centroid_rec2$distance) ## 25th = 22.73km, 75th = 67.32km

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
ggsave(rec_precip_centroid, file="Graphs/multidirectional_total_precip_rec_new.png")

## Correlate direction of actual range shift and temperature shift 
library(CircStats)
precip_radians <- precip_centroid_rec2[,c(1,13)]
colnames(precip_radians)[2] <- "precip_radians"
range_radians <- range_centroid_rec2[,c(1,13)]
colnames(range_radians)[2] <- "range_radians"
radians <- merge(precip_radians, range_radians, by="Common_name", all.x=TRUE)
circ.cor(radians$precip_radians, radians$range_radians, test=TRUE)
## significant - species are moving in a similar direction to what is expected for them
## to stay within their temperature envelope
## r = 0.9

## correlate the distances too
precip_dist <- precip_centroid_rec2[,c(1,10)]
colnames(precip_dist)[2] <- "precip_dist"
range_dist <- range_centroid_rec2[,c(1,10)]
colnames(range_dist)[2] <- "range_dist"
distance <- merge(precip_dist, range_dist, by="Common_name", all.x=TRUE)

plot(distance$precip_dist, distance$range_dist)
cor.test(distance$precip_dist, distance$range_dist, method="spearman")
## 0.94, p<0.001


## put all circular plots together
library(ggpubr)
circular_rec <- ggarrange(rec_centroid, rec_temp_centroid, rec_precip_centroid,
                          labels = c("(a)", "(b)", "(c)"),
                          ncol = 3, nrow = 1)
ggsave(circular_rec, file="Graphs/circlar_range_temp_precip_rec.png", height=12, width=10)

#### 3. Define colonised, persisted, and extirpated hectads
nmrsdata_rec$Recorded <- 1
nmrsdata_rec_expand <- nmrsdata_rec %>% expand(Common_name, Time_period, Hectad)
recorded <- nmrsdata_rec[,c("Common_name", "Time_period", "Hectad", "Recorded")]
nmrsdata_rec_expand <- merge(nmrsdata_rec_expand, recorded, by=c("Common_name", "Time_period", "Hectad"), all.x=TRUE)
nmrsdata_rec_expand[is.na(nmrsdata_rec_expand)] <- 0 ## change NAs to zero == species was NOT recorded at this hectad in this time period
## change to long format
nmrsdata_rec_expand2 <- nmrsdata_rec_expand %>%
  spread(Time_period, Recorded) ## each species has 1424 rows = the number of recorded hectads
## first remove rows where TP1 AND TP2 == 0 (this a hectad where a species was never recorded)
nmrsdata_rec_expand2<-nmrsdata_rec_expand2[!(nmrsdata_rec_expand2$`1975-1991`==0 & nmrsdata_rec_expand2$`2012-2016`==0),]
nmrsdata_rec_expand2$Hectad_category <- case_when(
  nmrsdata_rec_expand2$`1975-1991`==0 & nmrsdata_rec_expand2$`2012-2016`==1 ~ "Colonisation",
  nmrsdata_rec_expand2$`1975-1991`==1 & nmrsdata_rec_expand2$`2012-2016`==0 ~ "Extirpation",
  TRUE ~ "Persistence"
  )

## plot results on a map
## add in lat/lon info
lat_lon <- nmrsdata_rec[,c("Hectad", "lat", "lon","elevation10x10km","elevation10x10km_SD")]
lat_lon <- lat_lon %>% distinct(Hectad, .keep_all = TRUE) ## different way of unique
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
## remove colonisations - total species are those present in TP1 (extirpations and persistence's)
nmrsdata_rec_expand3<-nmrsdata_rec_expand2[!(nmrsdata_rec_expand2$Hectad_category=="Colonisation"),]
extirpated_hecs1 <- nmrsdata_rec_expand3 %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
                   dplyr::summarise(tot_sp = n()) 
extirpated_hecs2 <- nmrsdata_rec_expand3 %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
  filter(Hectad_category=="Extirpation") %>%
  dplyr::summarise(extir_sp = n())
extirpated_hecs <- merge(extirpated_hecs1, extirpated_hecs2, by=c("Hectad", "lat", "lon", "elevation10x10km", "elevation10x10km_SD"))
extirpated_hecs$extir_prop <- extirpated_hecs$extir_sp/extirpated_hecs$tot_sp ## 841 hectads which have been extirpated
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
ggsave(extir_hecs, file="Maps/Rec_extirpated_hecs_new.png")

## are hectads with lower proportion of extirpated species at higher elevations?
## prediction: species being lost from lower elevation as they are moving uphill
ggplot(extirpated_hecs, aes(x = elevation10x10km, y = extir_prop)) + 
  geom_point() +
  geom_smooth(color = "blue") +
  theme_bw() ## possible negative relationship

cor.test(extirpated_hecs$extir_prop, extirpated_hecs$elevation10x10km, method="spearman", exact=FALSE)

# #################################################################
# #### How does the climate and elevation of extirpated hectads compared to all recorded hectads?
# 
# ## proportion plots that Andy suggested
# ## subset 'well' extirpated hecs which have lost at least 25% of species
# well_extir_hecs <- extirpated_hecs[extirpated_hecs$extir_prop>=0.25,] ## 813 hectads
# well_extir <- well_extir_hecs$Hectad
# rec_hecs <- nmrsdata_rec[,c("Hectad","elevation10x10km","elevation10x10km_SD")] ## all rec hecs: 1424
# rec_hecs <- rec_hecs %>% distinct(Hectad, .keep_all = TRUE) ## different way of unique
# #rec_hecs <- rec_hecs %>% dplyr::filter(!Hectad %in% well_extir) ## now 611
# 
# ## add in climate data to each data frame
# temp <- read.csv("Data/NMRS/All_NMRS_hectads_annual_temperature.csv", header=TRUE)
# summer_temp <- read.csv("Data/NMRS/All_NMRS_hectads_summer_temperature.csv", header=TRUE)
# precip <- read.csv("Data/NMRS/All_NMRS_hectads_total_precipitation.csv", header=TRUE)
# ## put these together
# df_list <- list(temp, summer_temp, precip)
# climate <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)
# ## only want TP2 climate data
# climate_tp2 <- climate[climate$time_period=="TP2",]
# 
# ## merge in with all recorded hectads and then all well extirpated hectads
# rec_hecs <- merge(rec_hecs, climate_tp2, by="Hectad", all.x=)
# well_extir_hecs <- merge(well_extir_hecs, climate_tp2, by="Hectad", all.x=)
# 
# #### Annual temperature
# # split temperature into 10 bins
# rec_hecs$temperature_cat <- cut(rec_hecs$temperature, breaks = 10, labels=FALSE)
# 
# rec_hecs_temp_prop <- rec_hecs %>% group_by(temperature_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(rec_hecs_temp_prop$freq ~ rec_hecs_temp_prop$temperature_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(temperature_cat = cut(rec_hecs$temperature, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(temperature_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# min_max$cat <- 1:10
# 
# medians <- rec_hecs %>% group_by(temperature_cat) %>% dplyr::summarise(median=median(temperature))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(temperature_cat = case_when(
#   between(temperature, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(temperature, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(temperature, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(temperature, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(temperature, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(temperature, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(temperature, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(temperature, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(temperature, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(temperature_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$temperature_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("temperature_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# rec_hecs_temp_prop <- rec_hecs_temp_prop[,c("temperature_cat", "freq")]
# rec_hecs_temp_prop$hec_cat <- "All hectads"
# 
# temp_prop <- rbind(extirpated_hecs_prop, rec_hecs_temp_prop)
# temp_prop <- merge(temp_prop, medians, by="temperature_cat")
# ## plot line graph
# rec_extir_temp <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Mean annual temperature", y="Proportion of hectads") +
#   scale_x_continuous(limits=c(4,12)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# rec_extir_temp
# ggsave(rec_extir_temp, file="Graphs/Rec_hecs_extirpated_temperature.png")
# 
# 
# #### Summer temperature
# # split temperature into 10 bins
# rec_hecs$summer_temperature_cat <- cut(rec_hecs$summer_temperature, breaks = 10, labels=FALSE)
# 
# rec_hecs_temp_prop <- rec_hecs %>% group_by(summer_temperature_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(rec_hecs_temp_prop$freq ~ rec_hecs_temp_prop$summer_temperature_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(summer_temperature_cat = cut(rec_hecs$summer_temperature, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(summer_temperature_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# 
# medians <- rec_hecs %>% group_by(summer_temperature_cat) %>% dplyr::summarise(median=median(summer_temperature))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(summer_temperature_cat = case_when(
#   between(summer_temperature, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(summer_temperature, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(summer_temperature, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(summer_temperature, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(summer_temperature, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(summer_temperature, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(summer_temperature, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(summer_temperature, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(summer_temperature, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(summer_temperature_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$summer_temperature_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("summer_temperature_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# rec_hecs_temp_prop <- rec_hecs_temp_prop[,c("summer_temperature_cat", "freq")]
# rec_hecs_temp_prop$hec_cat <- "All hectads"
# 
# temp_prop <- rbind(extirpated_hecs_prop, rec_hecs_temp_prop)
# temp_prop <- merge(temp_prop, medians, by="summer_temperature_cat")
# ## plot line graph
# rec_extir_summer_temp <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Mean annual summer temperature", y="Proportion of hectads") +
#   scale_x_continuous(limits=c(10,18)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# rec_extir_summer_temp
# ggsave(rec_extir_summer_temp, file="Graphs/Rec_hecs_extirpated_summer_temperature.png")
# 
# 
# 
# #### Total precipitation
# # split precipitatino into 10 bins
# rec_hecs$precip_cat <- cut(rec_hecs$total_precip, breaks = 10, labels=FALSE)
# 
# rec_hecs_precip_prop <- rec_hecs %>% group_by(precip_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(rec_hecs_precip_prop$freq ~ rec_hecs_precip_prop$precip_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(precip_cat = cut(rec_hecs$total_precip, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(precip_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# 
# medians <- rec_hecs %>% group_by(precip_cat) %>% dplyr::summarise(median=median(total_precip))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(precip_cat = case_when(
#   between(total_precip, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(total_precip, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(total_precip, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(total_precip, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(total_precip, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(total_precip, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(total_precip, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(total_precip, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(total_precip, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(precip_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$precip_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("precip_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# rec_hecs_precip_prop <- rec_hecs_precip_prop[,c("precip_cat", "freq")]
# rec_hecs_precip_prop$hec_cat <- "All hectads"
# 
# precip_prop <- rbind(extirpated_hecs_prop, rec_hecs_precip_prop)
# precip_prop <- merge(precip_prop, medians, by="precip_cat")
# ## plot line graph
# rec_extir_precip <- ggplot(precip_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Mean total precipitation (mm)", y="Proportion of hectads") +
#   scale_x_continuous(breaks=seq(500,3500, by=500)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# rec_extir_precip
# ggsave(rec_extir_precip, file="Graphs/Rec_hecs_extirpated_total_precipitation.png")
# 
# 
# #### Elevation
# # split elevation into 10 bins
# rec_hecs$elev_cat <- cut(rec_hecs$elevation10x10km, breaks = 10, labels=FALSE)
# 
# rec_hecs_elev_prop <- rec_hecs %>% group_by(elev_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(rec_hecs_elev_prop$freq ~ rec_hecs_elev_prop$elev_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(elev_cat = cut(rec_hecs$elevation10x10km, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(elev_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# 
# medians <- rec_hecs %>% group_by(elev_cat) %>% dplyr::summarise(median=median(elevation10x10km))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(elev_cat = case_when(
#   between(elevation10x10km, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(elevation10x10km, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(elevation10x10km, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(elevation10x10km, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(elevation10x10km, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(elevation10x10km, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(elevation10x10km, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(elevation10x10km, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(elevation10x10km, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(elev_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt))
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$elev_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("elev_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# rec_hecs_elev_prop <- rec_hecs_elev_prop[,c("elev_cat", "freq")]
# rec_hecs_elev_prop$hec_cat <- "All hectads"
# 
# elev_prop <- rbind(extirpated_hecs_prop, rec_hecs_elev_prop)
# elev_prop <- merge(elev_prop, medians, by="elev_cat")
# ## plot line graph
# rec_extir_elev <- ggplot(elev_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Mean elevation (m)", y="Proportion of hectads") +
#   scale_x_continuous(breaks=seq(0,2000, by=200)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# rec_extir_elev
# ggsave(rec_extir_elev, file="Graphs/Rec_hecs_extirpated_elevation.png")
# 
# 
# #### SD Elevation
# # split elevation into 10 bins
# rec_hecs$elev_cat <- cut(rec_hecs$elevation10x10km_SD, breaks = 10, labels=FALSE)
# 
# rec_hecs_elev_prop <- rec_hecs %>% group_by(elev_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(rec_hecs_elev_prop$freq ~ rec_hecs_elev_prop$elev_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(elev_cat = cut(rec_hecs$elevation10x10km_SD, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(elev_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# 
# medians <- rec_hecs %>% group_by(elev_cat) %>% dplyr::summarise(median=median(elevation10x10km_SD))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(elev_cat = case_when(
#   between(elevation10x10km_SD, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(elevation10x10km_SD, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(elevation10x10km_SD, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(elevation10x10km_SD, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(elevation10x10km_SD, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(elevation10x10km_SD, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(elevation10x10km_SD, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(elevation10x10km_SD, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(elevation10x10km_SD, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(elev_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$elev_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("elev_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# rec_hecs_elev_prop <- rec_hecs_elev_prop[,c("elev_cat", "freq")]
# rec_hecs_elev_prop$hec_cat <- "All hectads"
# 
# elev_prop <- rbind(extirpated_hecs_prop, rec_hecs_elev_prop)
# elev_prop <- merge(elev_prop, medians, by="elev_cat")
# ## plot line graph
# rec_extir_elev_sd <- ggplot(elev_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Standard deviation elevation (m)", y="Proportion of hectads") +
#   scale_x_continuous(limits=c(0,450)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# rec_extir_elev_sd
# ggsave(rec_extir_elev_sd, file="Graphs/Rec_hecs_extirpated_elevation_SD.png")
# 
# 
# ### Are hectads which have warmed more over time more likely to be extirpated?
# ## calculate difference in climate (tp2 - tp1)
# temp_df <- climate[,c("Hectad", "time_period", "temperature")]
# climate_diff <- temp_df %>%
#   spread(time_period, temperature) %>% 
#   mutate(temp_diff = TP2-TP1)
# 
# ## merge in with all recorded hectads and then all well extirpated hectads
# rec_hecs <- merge(rec_hecs, climate_diff, by="Hectad", all.x=TRUE)
# well_extir_hecs <- merge(well_extir_hecs, climate_diff, by="Hectad", all.x=TRUE)
# 
# #### Annual temperature
# # split temperature into 10 bins
# rec_hecs$temp_diff_cat <- cut(rec_hecs$temp_diff, breaks = 10, labels=FALSE)
# 
# rec_hecs_temp_prop <- rec_hecs %>% group_by(temp_diff_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(rec_hecs_temp_prop$freq ~ rec_hecs_temp_prop$temp_diff_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(temp_diff_cat = cut(rec_hecs$temp_diff, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(temp_diff_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# 
# medians <- rec_hecs %>% group_by(temp_diff_cat) %>% dplyr::summarise(median=median(temp_diff))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(temp_diff_cat = case_when(
#   between(temp_diff, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(temp_diff, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(temp_diff, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(temp_diff, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(temp_diff, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(temp_diff, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(temp_diff, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(temp_diff, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(temp_diff, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(temp_diff_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$temp_diff_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("temp_diff_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# rec_hecs_temp_prop <- rec_hecs_temp_prop[,c("temp_diff_cat", "freq")]
# rec_hecs_temp_prop$hec_cat <- "All hectads"
# 
# temp_prop <- rbind(extirpated_hecs_prop, rec_hecs_temp_prop)
# temp_prop <- merge(temp_prop, medians, by="temp_diff_cat")
# ## plot line graph
# rec_extir_temp_diff <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Difference in temperature", y="Proportion of hectads") +
#   scale_x_continuous(limits=c(0.25,1.15)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# rec_extir_temp_diff
# ggsave(rec_extir_temp_diff, file="Graphs/Rec_hecs_extirpated_temperature_difference.png")
# 
# ## Put line graphs together
# climate_extir <- ggarrange(rec_extir_temp, rec_extir_summer_temp, rec_extir_temp_diff, rec_extir_precip,
#                            rec_extir_elev, rec_extir_elev_sd, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
#                          common.legend = TRUE, ncol = 2, nrow = 3, font.label = list(size = 12))
# ggsave(climate_extir, file="Graphs/Rec_climate_extirpated_hecs2.png", height=10, width=8)


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

# ## Are species that are moving southwards moving uphill?
# # Test by calculating mean elevation change between time periods
# mean_elev_well <- nmrsdata_well %>% group_by(Common_Name, time_period) %>%
#   dplyr::summarise(mean_elev=mean(elevation10x10km))
# ## change orientation
# mean_elev_well2 <- mean_elev_well %>%
#   gather(key, value, -Common_Name, -time_period) %>%
#   unite(col, key, time_period) %>%
#   spread(col, value)
# ## remove NAs - some species are only found in one TP 
# mean_elev_well2 <- na.omit(mean_elev_well2) ## 67 species
# 
# mean_elev_well2$elev_shift <- mean_elev_well2$mean_elev_TP2 - mean_elev_well2$mean_elev_TP1 ## TP2 - TP1
# hist(mean_elev_well2$elev_shift) ## looks good
# 
# mean(mean_elev_well2$elev_shift) ## 29.4 m uphill
# sd(mean_elev_well2$elev_shift) ## 142.5m 
# 
# # Shapiro-Wilk normality test for the differences
# shapiro.test(mean_elev_well2$elev_shift) # => p-value = 0.02 - distribution is NOT normal
# qqnorm(mean_elev_well2$elev_shift)
# qqline(mean_elev_well2$elev_shift) ## not good
# 
# ## Wilcox test to test whether elevation shifts differ significantly from 0
# wilcox.test(mean_elev_well2$elev_shift) ## significant
# 
# ## wilcoxon signed-rank test - does mean elevation differ between time periods?
# elev_well <- wilcox.test(mean_elev ~ time_period, data = mean_elev_well, paired = TRUE)
# elev_well ## non-significant
# 
# ## plot result
# elev_well_plot <- ggpaired(mean_elev_well, x = "time_period", y = "mean_elev",
#                           color = "time_period", line.color = "gray", line.size = 0.4,
#                           palette = "jco", id="Common_Name")+
#   xlab("Time period")+
#   ylab("Mean elevation")+
#   stat_compare_means(method="wilcox.test", paired = TRUE)
# elev_well_plot ## no change in mean elevation between TP1 and TP2
# ggsave(elev_well_plot, file="Graphs/Mean_elevation_shift_well_hecs.png")
# 
# length(which(mean_elev_well2$elev_shift>0)) ## 39 species moving uphill
# length(which(mean_elev_well2$elev_shift<0)) ## 26 species moving downhill
# ## categorise species as either moving uphill or downhill
# mean_elev_well2$shift_direction <- ifelse(mean_elev_well2$elev_shift>0, "Uphill", "Downhill")
# 
# ## add in elevation results to direction and distance shifts
# centroid_elev_well <- merge(mean_elev_well2, range_centroid_well2, by="Common_Name")
# 
# ## plot result coloured by elevational shift
# well_centroid2 <- ggplot(data=centroid_elev_well,
#                         aes(x=direction, y=distance, colour=shift_direction)) +
#   geom_segment(aes(xend = direction, yend = 0.1)) +
#   geom_point() +
#   scale_x_continuous(limits = c(0,360),
#                      breaks = seq(0, 360, by = 45),
#                      minor_breaks = seq(0, 360, by = 15)) +
#   # scale_y_log10() +
#   coord_polar() +
#   ggtitle("Well-recorded hectads") +
#   theme_bw()
# well_centroid2
# ggsave(well_centroid2, file="Graphs/multidirectional_elev_shifts_well.png")


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

# ###################
# # 2b: Mean summer temperature
# ###################
# 
# #### Direction and distance of shift based on mean annual temperature
# ## 2.5th and 97.5th percentiles for each species in TP1
# detach(package:Rmisc)
# detach(package:plyr)
# tp1_species_percentiles <- nmrsdata_well %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
#   summarise(summer_temp_2.5_perc = quantile(summer_temperature, probs = 0.025),
#             summer_temp_97.5_perc = quantile(summer_temperature, probs = 0.975))
# nmrsdata_well <- merge(nmrsdata_well, tp1_species_percentiles, by="Common_Name")
# 
# ## find hectads in TP2 which fall within 25th and 75th percentiles for each species
# filtered <- nmrsdata_well %>% group_by(Common_Name, time_period) %>%
#   filter(summer_temperature >= summer_temp_2.5_perc & summer_temperature <= summer_temp_97.5_perc)
# 
# ## calculate centroids of each
# summer_temp_centroid_well <- filtered %>% group_by(Common_Name, time_period) %>%
#   dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))
# 
# ## change orientation of data
# summer_temp_centroid_well2 <- summer_temp_centroid_well %>%
#   gather(key, value, -Common_Name, -time_period) %>%
#   unite(col, key, time_period) %>%
#   spread(col, value)
# ## NAs exist either because:
# # Species are only found in one hectad in TP1 => 25th & 75th percentile the same value, and no exact match in TP2
# # No hectads fall between 25th and 75th percentiles in TP2
# 
# ## remove NAs
# summer_temp_centroid_well2 <- na.omit(summer_temp_centroid_well2) ## 56 species
# 
# # then calculate distance (magnitude) and direction of shift between time periods for each species
# # Use pythagoras equation
# summer_temp_centroid_well2$distance <- sqrt((summer_temp_centroid_well2$mean_northing_TP1-summer_temp_centroid_well2$mean_northing_TP2)^2 + 
#                                              (summer_temp_centroid_well2$mean_easting_TP1-summer_temp_centroid_well2$mean_easting_TP2)^2) ## pythagoras equation
# summer_temp_centroid_well2$distance <- summer_temp_centroid_well2$distance/1000
# ## this gives distance in km
# 
# # use mean lat and lon values to calculate bearing/direction using geosphere package
# library(geosphere)
# summer_temp_centroid_well2$bearing <- bearing(summer_temp_centroid_well2[,c(6,4)], summer_temp_centroid_well2[,c(7,5)])
# summer_temp_centroid_well2$direction <- (summer_temp_centroid_well2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
# summer_temp_centroid_well2$radians <- NISTdegTOradian(summer_temp_centroid_well2$direction)
# x <- deg(circ.mean(summer_temp_centroid_well2$radians)) # -36, add full circle
# (x+360)%%360 ## 329
# ci <- vm.bootstrap.ci(summer_temp_centroid_well2$radians) 
# deg(c(5.45, 5.99)) ## lower=312 upper=343
# # median & quartiles for distances
# median(summer_temp_centroid_well2$distance) ## 43.7km 
# quantile(summer_temp_centroid_well2$distance) ## 25th = 25.1km, 75th = 75.8km
# 
# ## plot result
# ybreaks <- c(0,50,100,150,200)
# well_summer_temp_centroid <- ggplot(data=summer_temp_centroid_well2,
#                                    aes(x=direction, y=distance)) +
#   geom_segment(aes(xend = direction, yend = 0.1)) +
#   geom_point() +
#   scale_x_continuous(limits = c(0,360),
#                      breaks = seq(0, 360, by = 45),
#                      minor_breaks = seq(0, 360, by = 15)) +
#   coord_polar() +
#   scale_y_continuous(breaks = ybreaks) + 
#   geom_text(data = data.frame(x = 60, y = ybreaks, label = ybreaks),
#             aes(x = x, y = y, label = label),
#             inherit.aes = F,
#             size = 3) +
#   theme_bw() +
#   theme( panel.border = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
#          axis.title.x = element_blank(), axis.title.y = element_blank())
# well_summer_temp_centroid
# ggsave(well_summer_temp_centroid, file="Graphs/multidirectional_summer_temp_well.png")
# 
# ## Correlate direction of actual range shift and temperature shift 
# library(CircStats)
# temp_radians <- summer_temp_centroid_well2[,c("Common_Name", "radians")]
# colnames(temp_radians)[2] <- "temp_radians"
# range_radians <- range_centroid_well2[,c("Common_Name", "radians")]
# colnames(range_radians)[2] <- "range_radians"
# radians <- merge(temp_radians, range_radians, by="Common_Name", all.x=TRUE)
# circ.cor(radians$temp_radians, radians$range_radians, test=TRUE)
# ## significant - species are moving in a similar direction to what is expected for them
# ## to stay within their temperature envelope
# ## r = 0.68
# ## very similar to mean annual temperature

# 2c: Mean total precipitation
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

# 2c: Mean total precipitation
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






























# 
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# extirpated_hecs1 <- nmrsdata_well_expand2 %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
#   dplyr::summarise(tot_sp = n())
# extirpated_hecs2 <- nmrsdata_well_expand2 %>% group_by(Hectad, lat, lon, elevation10x10km, elevation10x10km_SD) %>%
#   filter(Hectad_category=="Extirpation") %>%
#   dplyr::summarise(extir_sp = n())
# extirpated_hecs2$group <- "Extirpated"
# extirpated_hecs <- merge(extirpated_hecs1, extirpated_hecs2, by=c("Hectad", "lat", "lon", "elevation10x10km", "elevation10x10km_SD"), all=TRUE)
# extirpated_hecs$extir_sp[is.na(extirpated_hecs$extir_sp)] <- 0
# extirpated_hecs$group[is.na(extirpated_hecs$group)] <- "Persisted"
# extirpated_hecs$extir_prop <- extirpated_hecs$extir_sp/extirpated_hecs$tot_sp ## 553
# ## very long way of doing this - gave up on finding a quicker way!
# 
# ## heatmap
# extir_hecs <- ggplot() +
#   geom_polygon(data = worldmap,
#                aes(x = long, y = lat, group = group),
#                fill = 'gray90', color = 'black') +
#   coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
#   theme_void() +
#   geom_point(data =extirpated_hecs,
#              aes(x = as.numeric(lon),
#                  y = as.numeric(lat), colour=extir_prop), size=1) +
#   scale_color_viridis_c(name="Proportion of \n extirpated species") + 
#   theme(title = element_text(size = 12))
# ## extirpated hectads with the most species lost tend to be in Northern Scotland/Pennines
# extir_hecs
# ggsave(extir_hecs, file="Maps/Well_extirpated_hecs.png")
# 
# 
# #################################################################
# #### How does the climate and elevation of extirpated hectads compared to all recorded hectads?
# 
# ## proportion plots that Andy suggested
# well_extir_hecs <- extirpated_hecs[extirpated_hecs$extir_prop>=0.5,] ## 472 hectads
# extir_hecs <- extirpated_hecs$Hectad
# well_hecs <- nmrsdata_well[,c("Hectad","elevation10x10km","elevation10x10km_SD")] ## all rec hecs: 639
# well_hecs <- well_hecs %>% distinct(Hectad, .keep_all = TRUE) ## different way of unique
# well_hecs <- well_hecs %>% dplyr::filter(!Hectad %in% extir_hecs) ## 167
# 
# 
# ## add in climate data to each data frame
# temp <- read.csv("Data/NMRS/All_NMRS_hectads_annual_temperature.csv", header=TRUE)
# summer_temp <- read.csv("Data/NMRS/All_NMRS_hectads_summer_temperature.csv", header=TRUE)
# precip <- read.csv("Data/NMRS/All_NMRS_hectads_total_precipitation.csv", header=TRUE)
# ## put these together
# df_list <- list(temp, summer_temp, precip)
# climate <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)
# ## only want TP2 climate data
# climate_tp2 <- climate[climate$time_period=="TP2",]
# 
# ## merge in with all recorded hectads and then all well extirpated hectads
# well_hecs <- merge(well_hecs, climate_tp2, by="Hectad", all.x=TRUE)
# well_extir_hecs <- merge(well_extir_hecs, climate_tp2, by="Hectad", all.x=TRUE)
# 
# extirpated_hecs <- merge(extirpated_hecs, climate_tp2, by="Hectad", all.x=TRUE)
# ## are hectads with lower proportion of extirpated species at higher elevations?
# ## prediction: species being lost from lower elevation as they are moving uphill
# temp <- ggplot(extirpated_hecs, aes(x = temperature, y = extir_prop)) + 
#   geom_point() +
#   geom_smooth(color = "black") +
#   theme_bw() ## possible positive relationship
# cor.test(extirpated_hecs$extir_prop, extirpated_hecs$temperature, method="spearman", exact=FALSE)
# 
# summer_temp <- ggplot(extirpated_hecs, aes(x = summer_temperature, y = extir_prop)) + 
#   geom_point() +
#   geom_smooth(color = "black") +
#   theme_bw() ## possible positive relationship
# 
# cor.test(extirpated_hecs$extir_prop, extirpated_hecs$summer_temperature, method="spearman", exact=FALSE)
# wilcox.test(extirpated_hecs$extir_prop, extirpated_hecs$summer_temperature)
# 
# library("ggpubr")
# 
# precip <- ggplot(extirpated_hecs, aes(x = total_precip, y = extir_prop)) + 
#   geom_point() +
#   geom_smooth(color = "black") +
#   theme_bw() ## possible negative relationship
# cor.test(extirpated_hecs$extir_prop, extirpated_hecs$total_precip, method="spearman", exact=FALSE)
# 
# temp_df <- climate[,c("Hectad", "time_period", "temperature")]
# climate_diff <- temp_df %>%
#   spread(time_period, temperature) %>% 
#   mutate(temp_diff = TP2-TP1)
# 
# ## merge in with all recorded hectads and then all well extirpated hectads
# extirpated_hecs <- merge(extirpated_hecs, climate_diff, by="Hectad", all.x=TRUE)
# extirpated_hecs <- subset(extirpated_hecs, select=-c(TP1,TP2,lat_centre,lon_centre,time_period))
# write.csv(extirpated_hecs, file="Data/Extirpated_hecs_climate_well.csv", row.names=FALSE)
# 
# temp_diff <- ggplot(extirpated_hecs, aes(x = temp_diff, y = extir_prop)) + 
#   geom_point() +
#   geom_smooth(color = "black") +
#   theme_bw() ## possible positive relationship
# cor.test(extirpated_hecs$extir_prop, extirpated_hecs$temp_diff, method="spearman", exact=FALSE)
# 
# elev <- ggplot(extirpated_hecs, aes(x = elevation10x10km, y = extir_prop)) + 
#   geom_point() +
#   geom_smooth(color = "black") +
#   theme_bw() ## possible positive relationship
# cor.test(extirpated_hecs$extir_prop, extirpated_hecs$elevation10x10km, method="spearman", exact=FALSE)
# 
# elev_sd <- ggplot(extirpated_hecs, aes(x = elevation10x10km_SD, y = extir_prop)) + 
#   geom_point() +
#   geom_smooth(color = "black") +
#   theme_bw() ## possible positive relationship
# cor.test(extirpated_hecs$extir_prop, extirpated_hecs$elevation10x10km_SD, method="spearman", exact=FALSE)
# 
# hist(extirpated_hecs$extir_prop) ## left skew
# extirpated_hecs$extir_prop_t <- (extirpated_hecs$extir_prop)^2 
# hist(extirpated_hecs$extir_prop_t) ## right skew
# 
# round(cor(extirpated_hecs[,c("temperature","summer_temperature","temp_diff","total_precip","elevation10x10km",
#                            "elevation10x10km_SD")]),3)
# ## 0.94 (annual temp ~ summer temp)
# ## -0.73 (annual temp ~ elevation)
# ## 0.7 (precipitation ~ SD elevation)
# ## 0.65 (summer temp ~ temp difference)
# ## -0.64 (summer temp ~ elevation)
# ## -0.62 (summer temp ~ precipitation)
# ## -0.6 (temp difference ~ elevation)
# 
# ## problematic variables: annual temperature
# 
# ### try some GLMMs without annual temp
# 
# # Which distribution?
# library(MASS)       
# 
# par(mfrow=c(2,2))
# car::qqp(extirpated_hecs$extir_prop, "norm")
# car::qqp(extirpated_hecs$extir_prop, "lnorm")
# nbinom <- fitdistr(extirpated_hecs$extir_prop, "Negative Binomial")
# car::qqp(extirpated_hecs$extir_prop, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
# poisson <- fitdistr(extirpated_hecs$extir_prop, "Poisson")
# car::qqp(extirpated_hecs$extir_prop, "pois", lambda=poisson$estimate)
# ## none look good
# 
# ## scale all variables
# extirpated_hecs$temperature_S <- scale(extirpated_hecs$temperature)
# extirpated_hecs$summer_temperature_S <- scale(extirpated_hecs$summer_temperature)
# extirpated_hecs$temp_diff_S <- scale(extirpated_hecs$temp_diff)
# extirpated_hecs$total_precip_S <- scale(extirpated_hecs$total_precip)
# 
# ## Try Poisson
# library(betareg)
# ## transform variable first (as we have ones in the dataset)
# ## if y also assumes the extremes 0 and 1, a useful transformation in practice is (y * (n−1) + 0.5) / n where n is the sample size.
# tot_sp = extirpated_hecs$tot_sp
# library(mgcv)
# model1 <- glm(extir_prop ~ scale(summer_temperature) + scale(temp_diff) + scale(total_precip) + scale(total_precip)*scale(temp_diff), data=extirpated_hecs,
#               family=binomial(link="logit"), weights=tot_sp, na.action = "na.fail")
# 
# par(mfrow=c(2,2))
# plot(model1)
# 
# summary(model1)
# par(mfrow=c(1,1))
# hist(residuals(model1))
# 
# plot(predict(model1),extirpated_hecs$extir_prop,
#      xlab="predicted",ylab="actual")
# abline(a=0,b=1)
# 
# 
# ## plot main effects 
# library(ggeffects)
# pred <- ggpredict(model1,se=TRUE,terms="summer_temperature")
# 
# ggplot(pred, aes(x, predicted)) +
#   geom_line() +
#   geom_point(data=extirpated_hecs, aes(summer_temperature, extir_prop)) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#   theme_classic()
# 
# pred <- ggpredict(model1,se=TRUE,terms="temp_diff")
# 
# ggplot(pred, aes(x, predicted)) +
#   geom_line() +
#   geom_point(data=extirpated_hecs, aes(temp_diff, extir_prop)) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#   theme_classic() ## non-significant
# 
# pred <- ggpredict(model1,se=TRUE,terms="total_precip")
# 
# ggplot(pred, aes(x, predicted)) +
#   geom_line() +
#   geom_point(data=extirpated_hecs, aes(total_precip, extir_prop)) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#   theme_classic()
# 
# pred <- ggpredict(model1,se=TRUE,terms="elevation10x10km")
# 
# ggplot(pred, aes(x, predicted)) +
#   geom_line() +
#   geom_point(data=extirpated_hecs, aes(elevation10x10km, extir_prop)) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
#   theme_classic()
# 
# ggpredict(model1, c("temp_diff", "total_precip")) %>% plot(rawdata=TRUE)
# pred <- ggpredict(model1, se=TRUE, c("temp_diff", "total_precip"))
# 
# 
# pred.dat = expand.grid(summer_temperature = mean(extirpated_hecs$summer_temperature),
#                        temp_diff = seq(min(extirpated_hecs$temp_diff), max(extirpated_hecs$temp_diff), length=100),
#                        total_precip = seq(min(extirpated_hecs$total_precip), max(extirpated_hecs$total_precip), length=100))
# 
# # Add predictions
# pred.dat$extir_prop = predict(model1, newdata=pred.dat, type="response")
# 
# ggplot(pred.dat, aes(temp_diff, total_precip, fill=extir_prop)) + 
#   geom_tile() +
#   scale_fill_gradient2(low="dodgerblue2", mid="white", high="red", 
#                        midpoint=mean(pred.dat$extir_prop)) +
#   geom_point(data=extirpated_hecs, aes(temp_diff, total_precip, colour=group)) +
#   labs(fill="Proportion of\n extirpated \nspecies") +
#   labs(x="Difference in annual temperature", y= "Total precipitation") +
#   theme_classic()
# #
# 
# 
# 
# library(MuMIn)
# d1 <- dredge(model1)
# 
# 
# 
# plots <- ggarrange(temp, summer_temp, temp_diff, precip, elev, elev_sd,
#                    ncol = 2, nrow = 3)
# ggsave(plots, file="Graphs/extir_prop_climate.png", height=12, width=10)
# #### Annual temperature
# # split temperature into 10 bins
# well_hecs$temperature_cat <- cut(well_hecs$temperature, breaks = 10, labels=FALSE)
# 
# well_hecs_temp_prop <- well_hecs %>% group_by(temperature_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt))
# plot(well_hecs_temp_prop$freq ~ well_hecs_temp_prop$temperature_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(temperature_cat = cut(well_hecs$temperature, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(temperature_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# min_max$cat <- 1:10
# 
# medians <- well_hecs %>% group_by(temperature_cat) %>% dplyr::summarise(median=median(temperature))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(temperature_cat = case_when(
#   between(temperature, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(temperature, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(temperature, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(temperature, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(temperature, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(temperature, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(temperature, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(temperature, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(temperature, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(temperature_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$temperature_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("temperature_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# well_hecs_temp_prop <- well_hecs_temp_prop[,c("temperature_cat", "freq")]
# well_hecs_temp_prop$hec_cat <- "All hectads"
# 
# temp_prop <- rbind(extirpated_hecs_prop, well_hecs_temp_prop)
# temp_prop <- merge(temp_prop, medians, by="temperature_cat")
# ## plot line graph
# well_extir_temp <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Mean annual temperature", y="Proportion of hectads") +
#   scale_x_continuous(limits=c(4,12)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# well_extir_temp
# ggsave(well_extir_temp, file="Graphs/Well_hecs_extirpated_temperature.png")
# 
# 
# #### Summer temperature
# # split temperature into 10 bins
# well_hecs$summer_temperature_cat <- cut(well_hecs$summer_temperature, breaks = 10, labels=FALSE)
# 
# well_hecs_temp_prop <- well_hecs %>% group_by(summer_temperature_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(well_hecs_temp_prop$freq ~ well_hecs_temp_prop$summer_temperature_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(summer_temperature_cat = cut(well_hecs$summer_temperature, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(summer_temperature_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# 
# medians <- well_hecs %>% group_by(summer_temperature_cat) %>% dplyr::summarise(median=median(summer_temperature))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(summer_temperature_cat = case_when(
#   between(summer_temperature, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(summer_temperature, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(summer_temperature, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(summer_temperature, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(summer_temperature, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(summer_temperature, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(summer_temperature, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(summer_temperature, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(summer_temperature, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(summer_temperature_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$summer_temperature_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("summer_temperature_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# well_hecs_temp_prop <- well_hecs_temp_prop[,c("summer_temperature_cat", "freq")]
# well_hecs_temp_prop$hec_cat <- "All hectads"
# 
# temp_prop <- rbind(extirpated_hecs_prop, well_hecs_temp_prop)
# temp_prop <- merge(temp_prop, medians, by="summer_temperature_cat")
# ## plot line graph
# well_extir_summer_temp <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Mean annual summer temperature", y="Proportion of hectads") +
#   scale_x_continuous(limits=c(10,18)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# well_extir_summer_temp
# ggsave(well_extir_summer_temp, file="Graphs/Well_hecs_extirpated_summer_temperature.png")
# 
# 
# 
# #### Total precipitation
# # split precipitatino into 10 bins
# well_hecs$precip_cat <- cut(well_hecs$total_precip, breaks = 10, labels=FALSE)
# 
# well_hecs_precip_prop <- well_hecs %>% group_by(precip_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(well_hecs_precip_prop$freq ~ well_hecs_precip_prop$precip_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(precip_cat = cut(well_hecs$total_precip, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(precip_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# 
# medians <- well_hecs %>% group_by(precip_cat) %>% dplyr::summarise(median=median(total_precip))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(precip_cat = case_when(
#   between(total_precip, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(total_precip, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(total_precip, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(total_precip, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(total_precip, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(total_precip, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(total_precip, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(total_precip, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(total_precip, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(precip_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$precip_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("precip_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# well_hecs_precip_prop <- well_hecs_precip_prop[,c("precip_cat", "freq")]
# well_hecs_precip_prop$hec_cat <- "All hectads"
# 
# precip_prop <- rbind(extirpated_hecs_prop, well_hecs_precip_prop)
# precip_prop <- merge(precip_prop, medians, by="precip_cat")
# ## plot line graph
# well_extir_precip <- ggplot(precip_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Mean total precipitation (mm)", y="Proportion of hectads") +
#   scale_x_continuous(breaks=seq(500,3500, by=500)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# well_extir_precip
# ggsave(well_extir_precip, file="Graphs/Well_hecs_extirpated_total_precipitation.png")
# 
# 
# #### Elevation
# # split elevation into 10 bins
# well_hecs$elev_cat <- cut(well_hecs$elevation10x10km, breaks = 10, labels=FALSE)
# 
# well_hecs_elev_prop <- well_hecs %>% group_by(elev_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(well_hecs_elev_prop$freq ~ well_hecs_elev_prop$elev_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(elev_cat = cut(well_hecs$elevation10x10km, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(elev_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# 
# medians <- well_hecs %>% group_by(elev_cat) %>% dplyr::summarise(median=median(elevation10x10km))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(elev_cat = case_when(
#   between(elevation10x10km, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(elevation10x10km, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(elevation10x10km, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(elevation10x10km, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(elevation10x10km, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(elevation10x10km, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(elevation10x10km, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(elevation10x10km, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(elevation10x10km, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(elev_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$elev_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("elev_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# well_hecs_elev_prop <- well_hecs_elev_prop[,c("elev_cat", "freq")]
# well_hecs_elev_prop$hec_cat <- "All hectads"
# 
# elev_prop <- rbind(extirpated_hecs_prop, well_hecs_elev_prop)
# elev_prop <- merge(elev_prop, medians, by="elev_cat")
# ## plot line graph
# well_extir_elev <- ggplot(elev_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Mean elevation (m)", y="Proportion of hectads") +
#   scale_x_continuous(breaks=seq(0,2000, by=200)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# well_extir_elev
# ggsave(well_extir_elev, file="Graphs/Well_hecs_extirpated_elevation.png")
# 
# 
# #### SD Elevation
# # split elevation into 10 bins
# well_hecs$elev_cat <- cut(well_hecs$elevation10x10km_SD, breaks = 10, labels=FALSE)
# 
# well_hecs_elev_prop <- well_hecs %>% group_by(elev_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(well_hecs_elev_prop$freq ~ well_hecs_elev_prop$elev_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(elev_cat = cut(well_hecs$elevation10x10km_SD, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(elev_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# 
# medians <- well_hecs %>% group_by(elev_cat) %>% dplyr::summarise(median=median(elevation10x10km_SD))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(elev_cat = case_when(
#   between(elevation10x10km_SD, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(elevation10x10km_SD, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(elevation10x10km_SD, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(elevation10x10km_SD, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(elevation10x10km_SD, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(elevation10x10km_SD, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(elevation10x10km_SD, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(elevation10x10km_SD, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(elevation10x10km_SD, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(elev_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$elev_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("elev_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# well_hecs_elev_prop <- well_hecs_elev_prop[,c("elev_cat", "freq")]
# well_hecs_elev_prop$hec_cat <- "All hectads"
# 
# elev_prop <- rbind(extirpated_hecs_prop, well_hecs_elev_prop)
# elev_prop <- merge(elev_prop, medians, by="elev_cat")
# ## plot line graph
# well_extir_elev_sd <- ggplot(elev_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Standard deviation elevation (m)", y="Proportion of hectads") +
#   scale_x_continuous(limits=c(0,400)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# well_extir_elev_sd
# ggsave(well_extir_elev_sd, file="Graphs/Well_hecs_extirpated_elevation_SD.png")
# 
# 
# ### Are hectads which have warmed more over time more likely to be extirpated?
# ## calculate difference in climate (tp2 - tp1)
# temp_df <- climate[,c("Hectad", "time_period", "temperature")]
# climate_diff <- temp_df %>%
#   spread(time_period, temperature) %>% 
#   mutate(temp_diff = TP2-TP1)
# 
# ## merge in with all recorded hectads and then all well extirpated hectads
# well_hecs <- merge(well_hecs, climate_diff, by="Hectad", all.x=TRUE)
# well_extir_hecs <- merge(well_extir_hecs, climate_diff, by="Hectad", all.x=TRUE)
# 
# #### Annual temperature
# # split temperature into 10 bins
# well_hecs$temp_diff_cat <- cut(well_hecs$temp_diff, breaks = 10, labels=FALSE)
# 
# well_hecs_temp_prop <- well_hecs %>% group_by(temp_diff_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(well_hecs_temp_prop$freq ~ well_hecs_temp_prop$temp_diff_cat) ## lower proportion of hectads with higher elevation
# ## same again for extirpated hectads
# 
# min_max <- tibble(temp_diff_cat = cut(well_hecs$temp_diff, breaks = 10)) %>% 
#   mutate(x_tmp = str_sub(temp_diff_cat, 2, -2)) %>% 
#   separate(x_tmp, c("min", "max"), sep = ",") %>% 
#   mutate_at(c("min", "max"), as.double)
# min_max <- unique(min_max)
# min_max <- arrange(min_max, min)
# 
# medians <- well_hecs %>% group_by(temp_diff_cat) %>% dplyr::summarise(median=median(temp_diff))
# 
# extirpated_hecs2 <- well_extir_hecs %>% mutate(temp_diff_cat = case_when(
#   between(temp_diff, min_max[1, "min"], min_max[1, "max"]) ~ 1,
#   between(temp_diff, min_max[2, "min"], min_max[2, "max"]) ~ 2,
#   between(temp_diff, min_max[3, "min"], min_max[3, "max"]) ~ 3,
#   between(temp_diff, min_max[4, "min"], min_max[4, "max"]) ~ 4,
#   between(temp_diff, min_max[5, "min"], min_max[5, "max"]) ~ 5,
#   between(temp_diff, min_max[6, "min"], min_max[6, "max"]) ~ 6,
#   between(temp_diff, min_max[7, "min"], min_max[7, "max"]) ~ 7,
#   between(temp_diff, min_max[8, "min"], min_max[8, "max"]) ~ 8,
#   between(temp_diff, min_max[9, "min"], min_max[9, "max"]) ~ 9,
#   TRUE ~ 10
# ))
# 
# extirpated_hecs_prop <- extirpated_hecs2 %>% group_by(temp_diff_cat) %>% 
#   dplyr::summarise(cnt = n()) %>%
#   mutate(freq = cnt / sum(cnt), 3)
# plot(extirpated_hecs_prop$freq ~ extirpated_hecs_prop$temp_diff_cat) ## lower proportion of hectads with higher elevation
# 
# extirpated_hecs_prop <- extirpated_hecs_prop[,c("temp_diff_cat", "freq")]
# extirpated_hecs_prop$hec_cat <- "Extirpated hectads"
# well_hecs_temp_prop <- well_hecs_temp_prop[,c("temp_diff_cat", "freq")]
# well_hecs_temp_prop$hec_cat <- "All hectads"
# 
# temp_prop <- rbind(extirpated_hecs_prop, well_hecs_temp_prop)
# temp_prop <- merge(temp_prop, medians, by="temp_diff_cat")
# ## plot line graph
# well_extir_temp_diff <- ggplot(temp_prop, aes(x=median, y=freq, group=hec_cat)) +
#   geom_line(aes(linetype=hec_cat), lwd=1) + 
#   labs(x="Difference in temperature", y="Proportion of hectads") +
#   scale_x_continuous(limits=c(0.4,1.1)) +
#   scale_y_continuous(breaks=seq(0,1, by=0.05)) +
#   theme_classic() +
#   theme(legend.title = element_blank())
# well_extir_temp_diff
# ggsave(well_extir_temp_diff, file="Graphs/Well_hecs_extirpated_temperature_difference.png")
# 
# 
# ## Put line graphs together
# climate_extir <- ggarrange(well_extir_temp, well_extir_summer_temp, well_extir_temp_diff, well_extir_precip,
#                            well_extir_elev, well_extir_elev_sd, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
#                            common.legend = TRUE, ncol = 2, nrow = 3, font.label = list(size = 12))
# ggsave(climate_extir, file="Graphs/Well_rec_climate_extirpated_hecs2.png", height=10, width=8)
# #
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
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

