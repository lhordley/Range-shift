##########################
#### user: Lisbeth Hordley
#### date: October 2021
#### info: Median & percentile temperature changes over time for upland moths

rm(list = ls())

## load packages
library(tidyverse)
library(lubridate)
library(ggplot2)

## NMRS temperature data for each TP
nmrsdata_temp_early <- readRDS("Data/NMRS/NMRS_hectad_temperature_TP1.rds") ## NMRS data for TP1 with temperature 
nmrsdata_temp_late <- readRDS("Data/NMRS/NMRS_hectad_temperature_TP2.rds") ## NMRS data for TP2 with temperature 
## migrant hectads that need excluded
migrant_hectads_early <- read.csv("Data/NMRS/migrant_hectad_exclusion_TP1.csv", header=TRUE) ## species x hectad combinations to remove due to immigrant populations
migrant_exclusion_late <- read.csv("Data/NMRS/migrant_hectad_exclusion_TP2.csv", header=TRUE)
## upland species list to filter
upland_species <- read.csv("Data/NMRS/NMRS_temp_median_percentile_selection.csv", header=TRUE)
## hectad recording levels
hec_records <- read.csv("Data/NMRS/Hectad_recording_levels_1975_1991_2012_2016.csv", header=TRUE)

### first exclude migrant hectads from each nmrs dataset
migrant_hectads_early[is.na(migrant_hectads_early)] <- 0
nmrsdata_temp_early <- merge(nmrsdata_temp_early, migrant_hectads_early, by=c("Common_Name", "Hectad"), all.x=TRUE)
nmrsdata_temp_early[is.na(nmrsdata_temp_early)] <- 0
nmrsdata_temp_early <- nmrsdata_temp_early[nmrsdata_temp_early$Exclude !=1, ] ## removes 84 rows (now 165390)

migrant_exclusion_late[is.na(migrant_exclusion_late)] <- 0
nmrsdata_temp_late <- merge(nmrsdata_temp_late, migrant_exclusion_late, by=c("Common_Name", "Hectad"), all.x=TRUE)
nmrsdata_temp_late[is.na(nmrsdata_temp_late)] <- 0
nmrsdata_temp_late <- nmrsdata_temp_late[nmrsdata_temp_late$Exclude !=1, ] ## removes 45 rows (now 390587)

nmrsdata_temp_early$time_period <- "TP1"
nmrsdata_temp_late$time_period <- "TP2"

## filter by upland/northern species only 
upland_species[is.na(upland_species)] <- 0 # replace NAs with 0 (0 = include the species as northern/upland)
upland_species_final <- upland_species[upland_species$exclude !=1, ] ## 73 species
upland <- upland_species_final[,1]
nmrsdata_temp_early <- nmrsdata_temp_early[which(nmrsdata_temp_early$Common_Name %in% upland), ] # 5390 rows
nmrsdata_temp_late <- nmrsdata_temp_late[which(nmrsdata_temp_late$Common_Name %in% upland), ] # 11079 rows
length(unique(nmrsdata_temp_early$Common_Name)) # 73 species
length(unique(nmrsdata_temp_late$Common_Name)) # 72 species
## which species is not in late time period? 
early_species <- data.frame(species=unique(nmrsdata_temp_early$Common_Name))
late_species <- data.frame(species=unique(nmrsdata_temp_late$Common_Name))
library(prodlim)
early_species$match <- ifelse(is.na(row.match(early_species, late_species)), "no", "yes")     
## scotch/mountain burnet is not present in TP2 (no records beyond 2011)
## remove from nmrsdata_temp_early
nmrsdata_temp_early <- nmrsdata_temp_early[!nmrsdata_temp_early$Common_Name=="Scotch Burnet or Mountain Burnet",]

## put datasets together
nmrsdata_temp_final <- rbind(nmrsdata_temp_early, nmrsdata_temp_late)


################ 
## split based on recording level of hectads

## Recorded hectads (hectads recorded in both time periods)
hectads <- unique(hec_records$HECTAD) ## 1782
nmrsdata_temp_rec <- nmrsdata_temp_final[which(nmrsdata_temp_final$Hectad %in% hectads), ]
length(unique(nmrsdata_temp_rec$Hectad)) ## 1424 hectads
length(unique(nmrsdata_temp_rec$Common_Name)) ## 72 species

detach(package:plyr)
library(dplyr)
summary_nmrs_temp <- nmrsdata_temp_early %>% 
  group_by(Common_Name) %>% 
  summarise(lower = min(temperature), upper = max(temperature), p = median(temperature))
summary_nmrs_temp <- summary_nmrs_temp[order(summary_nmrs_temp$p, decreasing = TRUE),]  

## test of differences in temperature between time periods
ggplot(nmrsdata_temp_rec, aes(x=time_period, y=temperature)) + 
  geom_boxplot() ## temp looks higher in TP2

library(Hmisc)
summary_nmrs_temp$groups<-as.numeric(cut2(summary_nmrs_temp$p, g=2))
## these groups are used to put groups into nmrs_temp to plot them more easily
sp_groups <- summary_nmrs_temp[,c(1,5)]
## use groups from above to put each species into groups to make plots
nmrsdata_temp_rec <- merge(nmrsdata_temp_rec, sp_groups, by="Common_Name", all=TRUE)
group <- unique(nmrsdata_temp_rec$groups) ## 2 groups
library(forcats)
nmrsdata_temp_rec$temperature <- as.numeric(as.character(nmrsdata_temp_rec$temperature))
legend_ord <- levels(with(nmrsdata_temp_rec, reorder(time_period, temperature)))

## plot map of sites for analysis: 1433 hectads and 72 upland species
worldmap = map_data('world')
rec_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_temp_final, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
rec_hecs
ggsave(rec_hecs, file="Maps/NMRS_cool_rec_hecs.png")

ggplot(data = nmrsdata_temp_rec[nmrsdata_temp_rec$groups==2 & nmrsdata_temp_rec$time_period=="TP1",], mapping = aes(x = reorder(Common_Name,-temperature, FUN=median), y = temperature)) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  #geom_hline(yintercept=9, linetype="dotted") +
  coord_flip() +
  labs(x="Common name", y="Temperature") +
  #scale_y_continuous(breaks = seq(3, 11, by = 1)) +
  scale_fill_discrete(breaks=legend_ord) +
  theme_light()
ggplot(data = nmrsdata_temp_rec[nmrsdata_temp_rec$groups==2,], mapping = aes(x = reorder(Common_Name,-temperature, FUN=median), y = temperature, fill=forcats::fct_rev(time_period))) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  #geom_hline(yintercept=9, linetype="dotted") +
  coord_flip() +
  labs(x="Common name", y="Temperature") +
  scale_fill_discrete(breaks=legend_ord) +
  #scale_y_continuous(breaks = seq(3, 11, by = 1)) +
  theme_light()
#

#### Direction and distance of shift based on mean annual temperature
## 25th and 75th percentiles for each species in TP1
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_temp_rec %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
  summarise(temp_25_perc = quantile(temperature, probs = 0.025),
            temp_75_perc = quantile(temperature, probs = 0.975))
nmrsdata_temp_rec <- merge(nmrsdata_temp_rec, tp1_species_percentiles, by="Common_Name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_temp_rec %>% group_by(Common_Name, time_period) %>%
  filter(temperature >= temp_25_perc & temperature <= temp_75_perc)

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
temp_centroid_rec2 <- na.omit(temp_centroid_rec2) ## 63 species (7)

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

## plot result
rec_centroid <- ggplot(data=temp_centroid_rec2,
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
ggsave(rec_centroid, file="Graphs/multidirectional_temp_rec.png")

## look at chestnut coloured carpet
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = nmrsdata_temp_rec[nmrsdata_temp_rec$Common_Name=="Angle-striped Sallow",], 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour=time_period), size=1) + 
  theme_void() +
  theme(title = element_text(size = 12))

### this method is very biased to recording effort - if the hectads with the temperature values
## between percentiles are not recorded, it won't give a true indication of what direction a 
## species should shift under climate change

## Well and heavily recorded hectads
well_heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded" | RECORDING.LEVEL=="Well recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
well_heavy_hecs <- well_heavy_hecs[which(well_heavy_hecs$n_row > 1), ]
well_heavy_hecs$n_row <- NULL ## 747 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_temp_wh <- nmrsdata_temp_final[which(nmrsdata_temp_final$Hectad %in% well_heavy_hecs$HECTAD), ]
length(unique(nmrsdata_temp_wh$Hectad)) ## 639 hectads
length(unique(nmrsdata_temp_wh$Common_Name)) ## 69 species

detach(package:plyr)
library(dplyr)
summary_nmrs_temp <- nmrsdata_temp_wh %>% 
  group_by(Common_Name) %>% 
  summarise(lower = min(temperature), upper = max(temperature), p = median(temperature)) 
summary_nmrs_temp <- summary_nmrs_temp[order(summary_nmrs_temp$p, decreasing = TRUE),]  
library(Hmisc)
summary_nmrs_temp$groups<-as.numeric(cut2(summary_nmrs_temp$p, g=2))
## these groups are used to put groups into nmrs_temp to plot them more easily
sp_groups <- summary_nmrs_temp[,c(1,5)] ## 69 species

# use groups from above to put each species into groups to make plots
nmrsdata_temp_wh <- merge(nmrsdata_temp_wh, sp_groups, by="Common_Name", all=TRUE)
group <- unique(nmrsdata_temp_wh$groups) ## 2 groups
library(forcats)
nmrsdata_temp_wh$temperature <- as.numeric(as.character(nmrsdata_temp_wh$temperature))
legend_ord <- levels(with(nmrsdata_temp_wh, reorder(time_period, temperature)))

## plot map of sites for analysis: 640 hectads and 69 upland species
worldmap = map_data('world')
wh_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_temp_wh, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
wh_hecs
ggsave(wh_hecs, file="Maps/NMRS_cool_well_hecs.png")

ggplot(data = nmrsdata_temp_wh[nmrsdata_temp_wh$groups==1,], mapping = aes(x = reorder(Common_Name,-temperature, FUN=median), y = temperature, fill=forcats::fct_rev(time_period))) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  #geom_hline(yintercept=9, linetype="dotted") +
  coord_flip() +
  labs(x="Common name", y="Temperature") +
  #scale_y_continuous(breaks = seq(3, 11, by = 1)) +
  scale_fill_discrete(breaks=legend_ord) +
  theme_light()
ggplot(data = nmrsdata_temp_wh[nmrsdata_temp_wh$groups==2,], mapping = aes(x = reorder(Common_Name,-temperature, FUN=median), y = temperature, fill=forcats::fct_rev(time_period))) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  #geom_hline(yintercept=9, linetype="dotted") +
  coord_flip() +
  labs(x="Common name", y="Temperature") +
  scale_fill_discrete(breaks=legend_ord) +
  #scale_y_continuous(breaks = seq(3, 11, by = 1)) +
  theme_light()
#

#### Direction and distance of shift based on mean annual temperature
## 25th and 75th percentiles for each species in TP1
detach(package:plyr)
tp1_species_percentiles <- nmrsdata_temp_wh %>% group_by(Common_Name) %>% filter(time_period=="TP1") %>%
  summarise(temp_25_perc = quantile(temperature, probs = 0.25),
            temp_75_perc = quantile(temperature, probs = 0.75))
nmrsdata_temp_wh <- merge(nmrsdata_temp_wh, tp1_species_percentiles, by="Common_Name")

## find hectads in TP2 which fall within 25th and 75th percentiles for each species
filtered <- nmrsdata_temp_wh %>% group_by(Common_Name, time_period) %>%
  filter(temperature >= temp_25_perc & temperature <= temp_75_perc)

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
temp_centroid_well2 <- na.omit(temp_centroid_well2) ## 50 species (lose 16 species)

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

## plot result
well_centroid <- ggplot(data=temp_centroid_well2,
                        aes(x=direction, y=distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  theme_bw()
well_centroid
ggsave(well_centroid, file="Graphs/multidirectional_temp_well.png")





## Heavily recorded hectads
heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
heavy_hecs <- heavy_hecs[which(heavy_hecs$n_row > 1), ]
heavy_hecs$n_row <- NULL ## 412 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_temp_heavy <- nmrsdata_temp_final[which(nmrsdata_temp_final$Hectad %in% heavy_hecs$HECTAD), ]
length(unique(nmrsdata_temp_heavy$Hectad)) ## 371 hectads

detach(package:plyr)
library(dplyr)
summary_nmrs_temp <- nmrsdata_temp_heavy %>% 
  group_by(Common_Name) %>% 
  summarise(lower = min(temperature), upper = max(temperature), p = median(temperature)) ## 796 species
summary_nmrs_temp <- summary_nmrs_temp[order(summary_nmrs_temp$p, decreasing = TRUE),]  
library(Hmisc)
summary_nmrs_temp$groups<-as.numeric(cut2(summary_nmrs_temp$p, g=2))
## these groups are used to put groups into nmrs_temp to plot them more easily
sp_groups <- summary_nmrs_temp[,c(1,5)] ## 64 species

# use groups from above to put each species into groups to make plots
nmrsdata_temp_heavy <- merge(nmrsdata_temp_heavy, sp_groups, by="Common_Name", all=TRUE)
group <- unique(nmrsdata_temp_heavy$groups) ## 2 groups
library(forcats)
nmrsdata_temp_heavy$temperature <- as.numeric(as.character(nmrsdata_temp_heavy$temperature))
legend_ord <- levels(with(nmrsdata_temp_heavy, reorder(time_period, temperature)))

## plot map of sites for analysis: 371 hectads and 64 upland species
worldmap = map_data('world')
heavy_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_temp_heavy[nmrsdata_temp_heavy$groups==1,], 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
heavy_hecs

ggplot(data = nmrsdata_temp_heavy[nmrsdata_temp_heavy$groups==1,], mapping = aes(x = reorder(Common_Name,-temperature, FUN=median), y = temperature, fill=forcats::fct_rev(time_period))) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  #geom_hline(yintercept=9, linetype="dotted") +
  coord_flip() +
  labs(x="Common name", y="Temperature") +
  #scale_y_continuous(breaks = seq(3, 11, by = 1)) +
  scale_fill_discrete(breaks=legend_ord) +
  theme_light()
ggplot(data = nmrsdata_temp_heavy[nmrsdata_temp_heavy$groups==2,], mapping = aes(x = reorder(Common_Name,-temperature, FUN=median), y = temperature, fill=forcats::fct_rev(time_period))) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  #geom_hline(yintercept=9, linetype="dotted") +
  coord_flip() +
  labs(x="Common name", y="Temperature") +
  scale_fill_discrete(breaks=legend_ord) +
  #scale_y_continuous(breaks = seq(3, 11, by = 1)) +
  theme_light()
#















