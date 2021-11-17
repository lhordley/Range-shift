##########################
#### user: Lisbeth Hordley
#### date: July 2021
#### info: Range margin shift analysis using NMRS data

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

### Filter: species need to occupy at least 20 hectads in each time period to calculate range shifts
# count the number of hectads per species for each year
nmrsdata_rec$PRESENCE <- 1
tp_species_hecs <- ddply(nmrsdata_rec, .(Preferred_Name,Common_Name, time_period), summarise,
                         HECTADS = sum(PRESENCE)) ## 144 (72 species in two time periods)
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_hecs <- tp_species_hecs %>% group_by(Common_Name) %>% filter(all(HECTADS>=20))
species_keep <- unique(tp_species_hecs$Common_Name) ## 44 species (lose 28)
## mostly removes rare/restricted species
nmrsdata_rec <-  nmrsdata_rec[which(nmrsdata_rec$Common_Name %in% species_keep),] ## 17 species
## remove duplicated records at 10km resolution - same species/year/hectad (probably have different 1km sites)
length(unique(nmrsdata_rec$Hectad)) ## 1422 hectads
length(unique(nmrsdata_rec$Common_Name)) ## 44 species

## map
worldmap = map_data('world')
rec_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_rec, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
rec_hecs
## save graph
#ggsave(rec_hecs, file="../../Maps/Upland_moth_17spp_rec_hecs.png")


################# LOW ELEVATION BOUNDARY SHIFT - mean elevation of 10 lowest elevation hectads for each time period

## order elevation by ascending value
## and take the top 10 values 
low_10_elev_rec_hecs <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>% arrange(elevation10x10km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_rec_hecs <- summarySE(low_10_elev_rec_hecs, measurevar="elevation10x10km", groupvars=c("Common_Name","time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_rec_hecs, 
          elevation10x10km[time_period == "TP1"] - elevation10x10km[time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = <0.001 - distribution is NOT normal
qqnorm(d)
qqline(d) ## not good

## wilcoxon signed-rank test instead
low_elev_rec <- wilcox.test(elevation10x10km ~ time_period, data = low_elev_rec_hecs, paired = TRUE)
low_elev_rec ## significant (p=0.007) - elevation is higher in TP1 compared to TP2

## plot result
rec_low_elev <- ggpaired(low_elev_rec_hecs, x = "time_period", y = "elevation10x10km",
                         color = "time_period", line.color = "gray", line.size = 0.4,
                         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
rec_low_elev
## low elevation boundary moving downhill over time
ggsave(rec_low_elev, file="Graphs/Low_elevation_shift_rec_hecs.png")

################# LOW LATITUDE BOUNDARY SHIFT - mean latitude of 10 most southerly hectads for each time period

## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_rec_hecs <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_rec_hecs <- summarySE(low_10_lat_rec_hecs, measurevar="lat", groupvars=c("Common_Name","time_period"))

## check for normality 
# compute the difference
d <- with(low_lat_rec_hecs, 
          lat[time_period == "TP1"] - lat[time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.6498 - normal distribution
qqnorm(d)
qqline(d) ## looks ok

## wilcoxon signed-rank test instead
low_lat_rec <- wilcox.test(lat ~ time_period, data = low_lat_rec_hecs, paired = TRUE)
low_lat_rec ## non-significant

## plot anyway
rec_low_lat <- ggpaired(low_lat_rec_hecs, x = "time_period", y = "lat",
                        color = "time_period", line.color = "gray", line.size = 0.4,
                        palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
rec_low_lat
ggsave(rec_low_lat, file="Graphs/Low_latitude_shift_rec_hecs.png")
## no change in low latitude boundary over time



##################################################################
##################### WELL RECORDED HECTADS ######################
##################################################################

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

########## Range margin shifts
####### Filter: species need to occupy at least 20 hectads in each time period to calculate range shifts
# count the number of hectads per species for each year
nmrsdata_well$PRESENCE <- 1
tp_species_hecs <- ddply(nmrsdata_well, .(Preferred_Name,Common_Name, time_period), summarise,
                         HECTADS = sum(PRESENCE)) ## 144 (72 species in two time periods)
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_hecs <- tp_species_hecs %>% group_by(Common_Name) %>% filter(all(HECTADS>=20))
species_keep <- unique(tp_species_hecs$Common_Name) ## 36 species
## mostly removes rare/restricted species
nmrsdata_well <-  nmrsdata_well[which(nmrsdata_well$Common_Name %in% species_keep),] ## 17 species
## remove duplicated records at 10km resolution - same species/year/hectad (probably have different 1km sites)
length(unique(nmrsdata_well$Hectad)) ## 639 hectads
length(unique(nmrsdata_well$Common_Name)) ## 36 species

## map
worldmap = map_data('world')
well_heavy_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_well, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
well_heavy_hecs
## save graph
#ggsave(well_heavy_hecs, file="../../Maps/Upland_moth_16spp_well_heavy_hecs.png")

################# LOW ELEVATION BOUNDARY SHIFT - mean elevation of 10 lowest elevation hectads for each time period

## order elevation by ascending value
## take the top 10 values 
low_10_elev_wh_hecs <- nmrsdata_well %>% group_by(Common_Name, time_period) %>% arrange(elevation10x10km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_wh_hecs <- summarySE(low_10_elev_wh_hecs, measurevar="elevation10x10km", groupvars=c("Common_Name","time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_wh_hecs, 
          elevation10x10km[time_period == "TP1"] - elevation10x10km[time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = <0.001 - distribution is NOT normal
qqnorm(d)
qqline(d) ## not good

## wilcoxon signed-rank test instead
low_elev_wh <- wilcox.test(elevation10x10km ~ time_period, data = low_elev_wh_hecs, paired = TRUE)
low_elev_wh ## non-significant

## plot anyway
low_elev_wh <- ggpaired(low_elev_wh_hecs, x = "time_period", y = "elevation10x10km",
                        color = "time_period", line.color = "gray", line.size = 0.4,
                        palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
low_elev_wh
ggsave(low_elev_wh, file="Graphs/Low_elevation_shift_well_hecs.png")


################# LOW LATITUDE BOUNDARY SHIFT - mean latitude of 10 most southerly hectads for each time period

## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_wh_hecs <- nmrsdata_well %>% group_by(Common_Name, time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_wh_hecs <- summarySE(low_10_lat_wh_hecs, measurevar="lat", groupvars=c("Common_Name","time_period"))

## check for normality 
# compute the difference
d <- with(low_lat_wh_hecs, 
          lat[time_period == "TP1"] - lat[time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.3095 - normal distribution
qqnorm(d)
qqline(d) ## not good

## wilcox.test
low_lat_wh <- wilcox.test(lat ~ time_period, data = low_lat_wh_hecs, paired = TRUE)
low_lat_wh ## non-significant

## plot anyway
low_lat_wh <- ggpaired(low_lat_wh_hecs, x = "time_period", y = "lat",
                       color = "time_period", line.color = "gray", line.size = 0.4,
                       palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="wilcox", paired = TRUE)
low_lat_wh
ggsave(low_lat_wh, file="Graphs/Low_latitude_shift_well_hecs.png")

