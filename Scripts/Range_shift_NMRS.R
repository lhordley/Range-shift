##########################
#### user: Lisbeth Hordley
#### date: July 2021
#### info: Range shift analysis using NMRS data

rm(list = ls())

## load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(sf)
library(ggplot2)
library(plyr)
library(ggpubr)
library(Rmisc)

## NMRS elevation data for each TP
nmrsdata <- readRDS("../../Range shift/Data/NMRS/NMRS_hectad_elevation.rds") ## NMRS data for all hectads and all years with elevation
## migrant hectads that need excluded
migrant_hectads_early <- read.csv("../../Range shift/Data/NMRS/migrant_hectad_exclusion_TP1.csv", header=TRUE) ## species x hectad combinations to remove due to immigrant populations
migrant_exclusion_late <- read.csv("../../Range shift/Data/NMRS/migrant_hectad_exclusion_TP2.csv", header=TRUE)
## upland species list to filter
upland_species <- read.csv("../../Range shift/Data/NMRS/NMRS_temp_median_percentile_selection.csv", header=TRUE)
## hectad recording levels
hec_records <- read.csv("../../Range shift/Data/NMRS/Hectad_recording_levels_1975_1991_2012_2016.csv", header=TRUE)
## temperature at all UK hectads in early & late TP
mean_temp_early <- read.csv("../../Range shift/Data/UKCP_tas_annual/Early_TP_75_91/mean_temp_75_91.csv", header=TRUE)
mean_temp_late <- read.csv("../../Range shift/Data/UKCP_tas_annual/Late_TP_12_16/mean_temp_12_16.csv", header=TRUE)

## subset nmrsdata by years of interest
early_years <- 1975:1991 
late_years <- 2012:2016
nmrsdata_early <- nmrsdata[which(nmrsdata$Year %in% early_years)]
nmrsdata_late <- nmrsdata[which(nmrsdata$Year %in% late_years)]
nmrsdata_early$time_period <- "TP1"
nmrsdata_late$time_period <- "TP2"

## take out migrant hectads
migrant_hectads_early[is.na(migrant_hectads_early)] <- 0
nmrsdata_early <- merge(nmrsdata_early, migrant_hectads_early, by=c("Common_Name", "Hectad"), all.x=TRUE)
nmrsdata_early[is.na(nmrsdata_early)] <- 0
nmrsdata_early <- nmrsdata_early[nmrsdata_early$Exclude !=1, ] 

migrant_exclusion_late[is.na(migrant_exclusion_late)] <- 0
nmrsdata_late <- merge(nmrsdata_late, migrant_exclusion_late, by=c("Common_Name", "Hectad"), all.x=TRUE)
nmrsdata_late[is.na(nmrsdata_late)] <- 0
nmrsdata_late <- nmrsdata_late[nmrsdata_late$Exclude !=1, ] 

## filter by upland/northern species only 
upland_species[is.na(upland_species)] <- 0 # replace NAs with 0 (0 = include the species as northern/upland)
upland_species_final <- upland_species[upland_species$exclude !=1, ] ## 73 species
upland <- upland_species_final[,1]
nmrsdata_early <- nmrsdata_early[which(nmrsdata_early$Common_Name %in% upland), ] 
nmrsdata_late <- nmrsdata_late[which(nmrsdata_late$Common_Name %in% upland), ] 
length(unique(nmrsdata_early$Common_Name)) # 73 species
length(unique(nmrsdata_late$Common_Name)) # 72 species
## which species is not in late time period? 
early_species <- data.frame(species=unique(nmrsdata_early$Common_Name))
late_species <- data.frame(species=unique(nmrsdata_late$Common_Name))
library(prodlim)
early_species$match <- ifelse(is.na(row.match(early_species, late_species)), "no", "yes")     
## scotch/mountain burnet is not present in TP2 (no records beyond 2011)
## remove from nmrsdata_temp_early
nmrsdata_early <- nmrsdata_early[!nmrsdata_early$Common_Name=="Scotch Burnet or Mountain Burnet",]

## put datasets together
nmrsdata <- rbind(nmrsdata_early, nmrsdata_late)
## remove year
nmrsdata$Year <- NULL
nmrsdata <- unique(nmrsdata) ## remove duplicated hectad/species combos within each time period
## 16373 rows

###############################################################
###################### RECORDED HECTADS #######################
###############################################################

## filter to recorded hectads only (i.e. those recorded once in both time periods)
hectads <- unique(hec_records$HECTAD) ## 1773
nmrsdata_rec <- nmrsdata[which(nmrsdata$Hectad %in% hectads), ]
length(unique(nmrsdata_rec$Hectad)) ## 1424 hectads
length(unique(nmrsdata_rec$Common_Name)) ## 72 species

# ## read in grey records - showing large southwards shift in range margin which is odd
# ## check raw nmrs records
# grey_recs <- read.csv("../../Range shift/Data/NMRS/nmrs_grey_records.csv", header=TRUE)
# grey_recs <- grey_recs[,c(5,8,14)]
# ## change gridref to 10km
# # remove white space
# grey_recs <- grey_recs %>% 
#   mutate(across(where(is.character), str_trim))
# grey_recs$gridref_10km <- dplyr::case_when(nchar(grey_recs$gridref)==10 ~ paste(substr(grey_recs$gridref,1,3), substr(grey_recs$gridref,7,7), sep=""), 
#                       nchar(grey_recs$gridref)==8 ~ paste(substr(grey_recs$gridref,1,3), substr(grey_recs$gridref,6,6), sep=""), 
#                       nchar(grey_recs$gridref)==6 ~ paste(substr(grey_recs$gridref,1,3), substr(grey_recs$gridref,5,5), sep=""),
#                                                                                               TRUE ~ grey_recs$gridref)
# 
# ## convert to lat/lon
# library(rnrfa)
# gridrefs <- as.data.frame(unique(grey_recs[,4])) ## 18 10km gridrefs
# colnames(gridrefs)[1] <- "Hectad"
# lon_lat <- as.data.frame(osg_parse(grid_refs=gridrefs$Hectad, coord_system = "WGS84"))
# gridref_latlon <- cbind(gridrefs, lon_lat)
# grey_recs <- merge(grey_recs, gridref_latlon, by.x="gridref_10km", by.y="Hectad")
# grey_recs$gridref <- NULL
# grey_recs <- unique(grey_recs) ## 70 unique 10km records
# ## plot grey moth
# years <- c(1975:1991, 2012:2016)
# worldmap = map_data('world')
# ggplot() + 
#   geom_polygon(data = worldmap, 
#                aes(x = long, y = lat, group = group), 
#                fill = 'gray90', color = 'black') + 
#   coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#   theme_void() + 
#   geom_point(data =grey_recs[grey_recs$rec_year %in% years,],
#              aes(x = as.numeric(lon),
#                  y = as.numeric(lat)), size=2) +
#   theme(title = element_text(size = 12))
# grey_recs <- grey_recs[grey_recs$rec_year %in% years,]
# grey <- nmrsdata[nmrsdata$Common_Name=="Grey",]
# worldmap = map_data('world')
# ggplot() + 
#   geom_polygon(data = worldmap, 
#                aes(x = long, y = lat, group = group), 
#                fill = 'gray90', color = 'black') + 
#   coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
#   theme_void() + 
#   geom_point(data =grey,
#              aes(x = as.numeric(lon),
#                  y = as.numeric(lat)), size=2) +
#   theme(title = element_text(size = 12))


##### Multidirectional centroid range shift
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
                       aes(direction, distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  ggtitle("Recorded hectads") +
  theme_bw()
rec_centroid

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

## has the temperature at species range centroids increase or decreased between time periods?
range_centroid_early <- range_centroid_rec[range_centroid_rec$time_period=="TP1",]
range_centroid_late <- range_centroid_rec[range_centroid_rec$time_period=="TP2",]

## now match temp locations with centroid locations in each time period

#### TP1
early_cent_latlon <- range_centroid_early[,c(3:4)]
early_cent_latlon <- data.frame(unique(early_cent_latlon)) ## 2732 rows

## need to source Find Closest Rcpp function first!
Rcpp::sourceCpp("Find_Closest_Rcpp_Function.cpp")

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


result <- find_closest(early_cent_latlon[, 1], early_cent_latlon[, 2], 
                       mean_temp_early[, 2], mean_temp_early[, 1])

result ### numbers are the corresponding row numbers that each centroid value matches in mean_temp_early?
## add row numbers to temp_latlon
early_cent_latlon$row_numbers <- rownames(early_cent_latlon)
early_cent_latlon$temp_row <- result
mean_temp_early$row_numbers <- rownames(mean_temp_early)

centroid_temp_tp1 <- merge(early_cent_latlon, mean_temp_early, by.x="temp_row", by.y="row_numbers")
centroid_temp_tp1 <- centroid_temp_tp1[,-c(1,4:6)]
range_centroid_early <- merge(range_centroid_early, centroid_temp_tp1, by=c("mean_lat", "mean_lon"))

#### TP1
late_cent_latlon <- range_centroid_late[,c(3:4)]
late_cent_latlon <- data.frame(unique(late_cent_latlon)) ## 2732 rows

result <- find_closest(late_cent_latlon[, 1], late_cent_latlon[, 2], 
                       mean_temp_late[, 2], mean_temp_late[, 1])

result ### numbers are the corresponding row numbers that each centroid value matches in mean_temp_early?
## add row numbers to temp_latlon
late_cent_latlon$row_numbers <- rownames(late_cent_latlon)
late_cent_latlon$temp_row <- result
mean_temp_late$row_numbers <- rownames(mean_temp_late)

centroid_temp_tp2 <- merge(late_cent_latlon, mean_temp_late, by.x="temp_row", by.y="row_numbers")
centroid_temp_tp2 <- centroid_temp_tp2[,-c(1,4:6)]
range_centroid_late <- merge(range_centroid_late, centroid_temp_tp2, by=c("mean_lat", "mean_lon"))

range_centroid_temp <- rbind(range_centroid_early, range_centroid_late)
## TP2 - TP1 temperature: positive values = temperature has increased between TP1 and TP2 and vice versa
range_centroid_temp2 <- range_centroid_temp %>% group_by(Common_Name) %>%
  dplyr::summarise(temp_centroid_diff = temperature[time_period == "TP2"] - temperature[time_period == "TP1"])
range_centroid_rec2 <- merge(range_centroid_rec2, range_centroid_temp2, by="Common_Name")

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

##### Range margin shifts
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
low_elev_rec ## significant (p=0.02) - elevation is higher in TP1 compared to TP2

## plot result
ggpaired(low_elev_rec_hecs, x = "time_period", y = "elevation10x10km",
         color = "time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
## low elevation boundary moving downhill over time

## method 2: accounting for status change using regression model
# calculate average latitude of 10 lowest grid cells in each time period = range margin
# range margin change = difference in range margin between time periods (TP2 - TP1)
# positive values == shift to north
## difference in low elevation boundary between TP1 and TP2
lrb_raw_shift <- low_elev_rec_hecs %>%
  group_by(Common_Name) %>%
  dplyr::summarise(lrb_diff = elevation10x10km[time_period == "TP2"] - elevation10x10km[time_period == "TP1"])
# negative values = srb has moved downhill
# positive values = srb has moved uphill (as predicted)

## ranges from 129m downhill (Light Knot Grass)
## to 461m uphill (Golden-rod Brindle)
## do these change after accounting for distribution change?


# distributional change = proportion of occupied hectads for each species out of all 
# hectads with observations of any species in TP1 and TP2
# change in distribution = TP2 - TP1
# no change in distribution = 0, positive = expansion, negative = contraction
all_hecs <- length(unique(nmrsdata_rec$Hectad))
status_change1 <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(prop_hectads = n()/all_hecs)
status_change <- status_change1 %>%
  group_by(Common_Name) %>%
  dplyr::summarise(status_diff = prop_hectads[time_period == "TP2"] - prop_hectads[time_period == "TP1"])
## beech-green carpet status_diff = 0 as no change in number of occupied hectads (70) between time periods

lrb_status <- merge(lrb_raw_shift, status_change, by="Common_Name")

# regression of change in change margin ~ change in distribution
# if intercept is positive and significantly different from zero, the species group has shifting SRB towards the north
# intercept is value of SRB when status_diff is 0
hist(lrb_status$lrb_diff) ## looks ok - one outlier (probably golden rod brindle) - might need removed (check cook's distance)
## looks good once golden-rod brindle is removed
lrb_status_model <- lm(lrb_diff ~ status_diff, data=lrb_status)
summary(lrb_status_model) 
## both non-significant - still non-significant after removing golden-rod brindle

## remove golden rod brindle and run model again
lrb_status <- lrb_status[!lrb_status$Common_Name=="Golden-rod Brindle",]

ggplot(lrb_status, aes(x=status_diff, y=lrb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in prevalence", y="Change in low elevation range margin (m)") +
  theme_classic()

# extract residuals to obtain corrected measure of SRB shifts per species 
# this is the residual shift in SRB that is not explain by the linear effect of prevalence across the studied species
# probably correlate with raw SRB estimates, but are more conservative (can show this in plot)
lrb_status_model$residuals        
# add these to srb_status dataframe
lrb_status$lrb_diff_adjust <- lrb_status_model$residuals
## plot relationship between the two
ggplot(lrb_status, aes(x=lrb_diff_adjust, y=lrb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in LRB (residuals)", y="Change in LRB (raw)") +
  theme_classic()
## very similar

## could split species into categories to show direction of responses (as per Hallfors et al 2021)
## then could do chi-squared between groups? Do species moving north also move uphill?

# srb_status$srb_category <- case_when(
#   srb_status$srb_diff_adjust>=20 ~ "Northwards",
#   srb_status$srb_diff_adjust<=-20 ~ "Southwards",
#   TRUE ~ "Stable"
# )
# proportions <- table(srb_status$srb_category)/length(srb_status$srb_category)
# percentages <- proportions*100

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
ggpaired(low_lat_rec_hecs, x = "time_period", y = "lat",
         color = "time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
## no change in low latitude boundary over time


## method 2: accounting for status change using regression model
# calculate average latitude of 10 southern-most grid cells in each time period = range margin
# range margin change = difference in range margin between time periods (TP2 - TP1)
# positive values == shift to north
low_10_northing_rec_hecs <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>% arrange(northing) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_northing_rec_hecs <- summarySE(low_10_northing_rec_hecs, measurevar="northing", groupvars=c("Common_Name","time_period"))
## difference in low latitude boundary between TP1 and TP2
srb_raw_shift <- low_northing_rec_hecs %>%
  group_by(Common_Name) %>%
  dplyr::summarise(srb_diff = northing[time_period == "TP2"] - northing[time_period == "TP1"])
# negative values = srb has moved southwards
# positive values = srb has moved northwards (as predicted)
## divide each by 1000 to give values in km - then residuals from model will also be in km
srb_raw_shift$srb_diff <- srb_raw_shift$srb_diff/1000
## ranges from 169km northwards (Pretty Pinion)
## to 105km southwards (Smoky Wave)
## do these change after accounting for distribution change?


# distributional change = proportion of occupied hectads for each species out of all 
# hectads with observations of any species in TP1 and TP2
# change in distribution = TP2 - TP1
# no change in distribution = 0, positive = expansion, negative = contraction
all_hecs <- length(unique(nmrsdata_rec$Hectad))
status_change1 <- nmrsdata_rec %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(prop_hectads = n()/all_hecs)
status_change <- status_change1 %>%
  group_by(Common_Name) %>%
  dplyr::summarise(status_diff = prop_hectads[time_period == "TP2"] - prop_hectads[time_period == "TP1"])
## beech-green carpet status_diff = 0 as no change in number of occupied hectads (70) between time periods

srb_status <- merge(srb_raw_shift, status_change, by="Common_Name")

# regression of change in change margin ~ change in distribution
# if intercept is positive and significantly different from zero, the species group has shifting SRB towards the north
# intercept is value of SRB when status_diff is 0
hist(srb_status$srb_diff) ## looks good
srb_status_model <- lm(srb_diff ~ status_diff, data=srb_status)
summary(srb_status_model) 
## intercept is significant & positive = species have overall shift SRB northwards than expected purely from their change in prevalence
## when x=0, the mean shift in SRB for species with no overall change in prevalence is a northward shift of 21.92km (p=0.038)
## shift in prevalence just non-significant (p=0.06) & negative

ggplot(srb_status, aes(x=status_diff, y=srb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in prevalence", y="Change in southern range margin (km)") +
  theme_classic()

# extract residuals to obtain corrected measure of SRB shifts per species 
# this is the residual shift in SRB that is not explain by the linear effect of prevalence across the studied species
# probably correlate with raw SRB estimates, but are more conservative (can show this in plot)
srb_status_model$residuals        
# add these to srb_status dataframe
srb_status$srb_diff_adjust <- srb_status_model$residuals
## plot relationship between the two
ggplot(srb_status, aes(x=srb_diff_adjust, y=srb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in SRB (residuals)", y="Change in SRB (raw)") +
  theme_classic()
## very similar

## can split species into categories to show direction of responses (as per Hallfors et al 2021)
# Northwards: shift of >= 20km
# Stable: shift beteween 20km and -20km
# Southwards: shift of <= -20km

srb_status$srb_category <- case_when(
  srb_status$srb_diff_adjust>=20 ~ "Northwards",
  srb_status$srb_diff_adjust<=-20 ~ "Southwards",
  TRUE ~ "Stable"
)
proportions <- table(srb_status$srb_category)/length(srb_status$srb_category)
percentages <- proportions*100


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

## plot result
wh_centroid <- ggplot(data=range_centroid_wh2,
                      aes(direction, distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  ggtitle("Well/heavily recorded hectads") +
  theme_bw()
wh_centroid

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

########## Range margin shifts
####### Filter: species need to occupy at least 20 hectads in each time period to calculate range shifts
# count the number of hectads per species for each year
nmrsdata_wh$PRESENCE <- 1
tp_species_hecs <- ddply(nmrsdata_wh, .(Preferred_Name,Common_Name, time_period), summarise,
                         HECTADS = sum(PRESENCE)) ## 144 (72 species in two time periods)
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_hecs <- tp_species_hecs %>% group_by(Common_Name) %>% filter(all(HECTADS>=20))
species_keep <- unique(tp_species_hecs$Common_Name) ## 36 species
## mostly removes rare/restricted species
nmrsdata_wh <-  nmrsdata_wh[which(nmrsdata_wh$Common_Name %in% species_keep),] ## 17 species
## remove duplicated records at 10km resolution - same species/year/hectad (probably have different 1km sites)
length(unique(nmrsdata_wh$Hectad)) ## 639 hectads
length(unique(nmrsdata_wh$Common_Name)) ## 36 species

## map
worldmap = map_data('world')
well_heavy_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_wh, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
well_heavy_hecs
## save graph
#ggsave(well_heavy_hecs, file="../../Maps/Upland_moth_16spp_well_heavy_hecs.png")

################# LOW ELEVATION BOUNDARY SHIFT - mean elevation of 10 lowest elevation hectads for each time period

## order elevation by ascending value
## take the top 10 values 
low_10_elev_wh_hecs <- nmrsdata_wh %>% group_by(Common_Name, time_period) %>% arrange(elevation10x10km) %>% slice(seq_len(10))
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
ggpaired(low_elev_wh_hecs, x = "time_period", y = "elevation10x10km",
         color = "time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)

## method 2: accounting for status change using regression model
# calculate average latitude of 10 lowest grid cells in each time period = range margin
# range margin change = difference in range margin between time periods (TP2 - TP1)
# positive values == shift to north
## difference in low elevation boundary between TP1 and TP2
lrb_raw_shift <- low_elev_wh_hecs %>%
  group_by(Common_Name) %>%
  dplyr::summarise(lrb_diff = elevation10x10km[time_period == "TP2"] - elevation10x10km[time_period == "TP1"])
# negative values = srb has moved downhill
# positive values = srb has moved uphill (as predicted)

## ranges from 129m downhill (Light Knot Grass)
## to 461m uphill (Golden-rod Brindle)
## do these change after accounting for distribution change?


# distributional change = proportion of occupied hectads for each species out of all 
# hectads with observations of any species in TP1 and TP2
# change in distribution = TP2 - TP1
# no change in distribution = 0, positive = expansion, negative = contraction
all_hecs <- length(unique(nmrsdata_wh$Hectad))
status_change1 <- nmrsdata_wh %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(prop_hectads = n()/all_hecs)
status_change <- status_change1 %>%
  group_by(Common_Name) %>%
  dplyr::summarise(status_diff = prop_hectads[time_period == "TP2"] - prop_hectads[time_period == "TP1"])
## beech-green carpet status_diff = 0 as no change in number of occupied hectads (70) between time periods

lrb_status <- merge(lrb_raw_shift, status_change, by="Common_Name")

# regression of change in change margin ~ change in distribution
# if intercept is positive and significantly different from zero, the species group has shifting SRB towards the north
# intercept is value of SRB when status_diff is 0
hist(lrb_status$lrb_diff) ## looks ok - few outliers - might need removed (check cook's distance)
## looks good once golden-rod brindle is removed
lrb_status_model <- lm(lrb_diff ~ status_diff, data=lrb_status)
summary(lrb_status_model) 
## intercept non-significant
## status_diff significant negative

ggplot(lrb_status, aes(x=status_diff, y=lrb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in prevalence", y="Change in low elevation range margin (m)") +
  theme_classic()

# extract residuals to obtain corrected measure of SRB shifts per species 
# this is the residual shift in SRB that is not explain by the linear effect of prevalence across the studied species
# probably correlate with raw SRB estimates, but are more conservative (can show this in plot)
lrb_status_model$residuals        
# add these to srb_status dataframe
lrb_status$lrb_diff_adjust <- lrb_status_model$residuals
## plot relationship between the two
ggplot(lrb_status, aes(x=lrb_diff_adjust, y=lrb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in LRB (residuals)", y="Change in LRB (raw)") +
  theme_classic()
## very similar

## could split species into categories to show direction of responses (as per Hallfors et al 2021)
## then could do chi-squared between groups? Do species moving north also move uphill?

# srb_status$srb_category <- case_when(
#   srb_status$srb_diff_adjust>=20 ~ "Northwards",
#   srb_status$srb_diff_adjust<=-20 ~ "Southwards",
#   TRUE ~ "Stable"
# )
# proportions <- table(srb_status$srb_category)/length(srb_status$srb_category)
# percentages <- proportions*100

################# LOW LATITUDE BOUNDARY SHIFT - mean latitude of 10 most southerly hectads for each time period
## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_wh_hecs <- nmrsdata_wh %>% group_by(Common_Name, time_period) %>% arrange(lat) %>% slice(seq_len(10))
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
ggpaired(low_lat_wh_hecs, x = "time_period", y = "lat",
         color = "time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="t.test", paired = TRUE)


## method 2: accounting for status change using regression model
# calculate average latitude of 10 southern-most grid cells in each time period = range margin
# range margin change = difference in range margin between time periods (TP2 - TP1)
# positive values == shift to north
low_10_northing_wh_hecs <- nmrsdata_wh %>% group_by(Common_Name, time_period) %>% arrange(northing) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_northing_wh_hecs <- summarySE(low_10_northing_wh_hecs, measurevar="northing", groupvars=c("Common_Name","time_period"))
## difference in low latitude boundary between TP1 and TP2
srb_raw_shift <- low_northing_wh_hecs %>%
  group_by(Common_Name) %>%
  dplyr::summarise(srb_diff = northing[time_period == "TP2"] - northing[time_period == "TP1"])
# negative values = srb has moved southwards
# positive values = srb has moved northwards (as predicted)
## divide each by 1000 to give values in km - then residuals from model will also be in km
srb_raw_shift$srb_diff <- srb_raw_shift$srb_diff/1000
## ranges from 169km northwards (Pretty Pinion)
## to 105km southwards (Smoky Wave)
## do these change after accounting for distribution change?


# distributional change = proportion of occupied hectads for each species out of all 
# hectads with observations of any species in TP1 and TP2
# change in distribution = TP2 - TP1
# no change in distribution = 0, positive = expansion, negative = contraction
all_hecs <- length(unique(nmrsdata_wh$Hectad))
status_change1 <- nmrsdata_wh %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(prop_hectads = n()/all_hecs)
status_change <- status_change1 %>%
  group_by(Common_Name) %>%
  dplyr::summarise(status_diff = prop_hectads[time_period == "TP2"] - prop_hectads[time_period == "TP1"])

srb_status <- merge(srb_raw_shift, status_change, by="Common_Name")

# regression of change in change margin ~ change in distribution
# if intercept is positive and significantly different from zero, the species group has shifting SRB towards the north
# intercept is value of SRB when status_diff is 0
hist(srb_status$srb_diff) ## looks good
srb_status_model <- lm(srb_diff ~ status_diff, data=srb_status)
summary(srb_status_model) 
## intercept is non-significant (16km, p=0.08)
## shift in prevalence just significant (p=0.008) & negative

ggplot(srb_status, aes(x=status_diff, y=srb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in prevalence", y="Change in southern range margin (km)") +
  theme_classic()

# extract residuals to obtain corrected measure of SRB shifts per species 
# this is the residual shift in SRB that is not explain by the linear effect of prevalence across the studied species
# probably correlate with raw SRB estimates, but are more conservative (can show this in plot)
srb_status_model$residuals        
# add these to srb_status dataframe
srb_status$srb_diff_adjust <- srb_status_model$residuals
## plot relationship between the two
ggplot(srb_status, aes(x=srb_diff_adjust, y=srb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in SRB (residuals)", y="Change in SRB (raw)") +
  theme_classic()
## very similar

## can split species into categories to show direction of responses (as per Hallfors et al 2021)
# Northwards: shift of >= 20km
# Stable: shift beteween 20km and -20km
# Southwards: shift of <= -20km

srb_status$srb_category <- case_when(
  srb_status$srb_diff_adjust>=20 ~ "Northwards",
  srb_status$srb_diff_adjust<=-20 ~ "Southwards",
  TRUE ~ "Stable"
)
proportions <- table(srb_status$srb_category)/length(srb_status$srb_category)
percentages <- proportions*100


###############################################################
## heavily recorded hectads
## filter to recorded hectads only (i.e. those recorded once in both time periods)
heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
heavy_hecs <- heavy_hecs[which(heavy_hecs$n_row > 1), ]
heavy_hecs$n_row <- NULL ## 411 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_heavy <- nmrsdata[which(nmrsdata$Hectad %in% heavy_hecs$HECTAD), ]
length(unique(nmrsdata_heavy$Hectad)) ## 370 hectads
length(unique(nmrsdata_heavy$Common_Name)) ## 64 species

#### Multidirection centroid range shift
# take mean easting and northing and lat and lon across each species occupied hectads for each time period = range centroid
# following Gillings et al 2015
range_centroid_heavy <- nmrsdata_heavy %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(mean_lat=mean(lat), mean_lon=mean(lon), mean_easting=mean(easting), mean_northing=mean(northing))
## change orientation of data
range_centroid_heavy2 <- range_centroid_heavy %>%
  gather(key, value, -Common_Name, -time_period) %>%
  unite(col, key, time_period) %>%
  spread(col, value)
## remove NAs - some species are only found in one TP 
range_centroid_heavy2 <- na.omit(range_centroid_heavy2) ## 59 species
# then calculate distance (magnitude) and direction of shift between time periods for each species
# Use pythagoras equation
range_centroid_heavy2$distance <- sqrt((range_centroid_heavy2$mean_northing_TP1-range_centroid_heavy2$mean_northing_TP2)^2 + 
                                         (range_centroid_heavy2$mean_easting_TP1-range_centroid_heavy2$mean_easting_TP2)^2) ## pythagoras equation
range_centroid_heavy2$distance <- range_centroid_heavy2$distance/1000
## this gives distance in km

# use mean lat and lon values to calculate bearing/direction using geosphere package
library(geosphere)
range_centroid_heavy2$bearing <- bearing(range_centroid_heavy2[,c(6,4)], range_centroid_heavy2[,c(7,5)])
range_centroid_heavy2$direction <- (range_centroid_heavy2$bearing + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360

## plot result
heavy_centroid <- ggplot(data=range_centroid_heavy2,
                         aes(direction, distance)) +
  geom_segment(aes(xend = direction, yend = 0.1)) +
  geom_point() +
  scale_x_continuous(limits = c(0,360),
                     breaks = seq(0, 360, by = 45),
                     minor_breaks = seq(0, 360, by = 15)) +
  # scale_y_log10() +
  coord_polar() +
  ggtitle("Heavily recorded hectads") +
  theme_bw()
heavy_centroid

## put plots together
library(gridExtra)
range_centroid <- grid.arrange(rec_centroid, wh_centroid, heavy_centroid, nrow = 1)
ggsave("../../Range shift/Graphs/NMRS_centroid_shift.png", range_centroid, width=15, height=10)

## check yellow-ringed carpet - moving far southwards
worldmap = map_data('world')
ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  theme_void() +
  geom_point(data =nmrsdata_rec[nmrsdata_rec$Common_Name=="Yellow-ringed Carpet",],
             aes(x = lon, y=lat, colour=time_period), size=1) +
  theme(title = element_text(size = 12))
## should species only be included if they occur in at least 2 hectads in each TP?
## or maybe can't trust heavily recorded hectads only?
x <- nmrsdata[nmrsdata$Common_Name=="Yellow-ringed Carpet",]

# use circ.mean from the circstats package to find a mean direction & confidence limits across all species
library(NISTunits)
library(CircStats)
range_centroid_heavy2$radians <- NISTdegTOradian(range_centroid_heavy2$direction)
x <- deg(circ.mean(range_centroid_heavy2$radians)) # -42, add full circle
(x+360)%%360 ## 314
# median & quartiles for distances
median(range_centroid_heavy2$distance) ## 37.7km 
quantile(range_centroid_heavy2$distance) ## 25th = 26.9km, 75th = 49.7km

## Rayleigh's test
## if significant, the distribution of directions differ significantly from an even distribution
## if R = 0 (the data are completely spread around the circle)
## if R = 1 (the data are completely concentrated on one point)
r <- r.test(range_centroid_heavy2$radians)
## p = 0.27
## r = 0.15
z <- (nrow(range_centroid_heavy2))*(r$r.bar^2) ## 1.29

data <- rvm(25, pi, 2)
r.test(data)

#### Range margin shifts
####### Filter: species need to occupy at least 20 hectads in each time period to calculate range shifts
# count the number of hectads per species for each year
nmrsdata_heavy$PRESENCE <- 1
tp_species_hecs <- ddply(nmrsdata_heavy, .(Preferred_Name,Common_Name, time_period), summarise,
                         HECTADS = sum(PRESENCE)) ## 144 (72 species in two time periods)
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_hecs <- tp_species_hecs %>% group_by(Common_Name) %>% filter(all(HECTADS>=20))
species_keep <- unique(tp_species_hecs$Common_Name) ## 29 species
## mostly removes rare/restricted species
nmrsdata_heavy <-  nmrsdata_heavy[which(nmrsdata_heavy$Common_Name %in% species_keep),] 
## remove duplicated records at 10km resolution - same species/year/hectad (probably have different 1km sites)
length(unique(nmrsdata_heavy$Hectad)) ## 367 hectads
length(unique(nmrsdata_heavy$Common_Name)) ## 29 species

## map
worldmap = map_data('world')
heavy_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_heavy, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
heavy_hecs
## save graph
## scottish coverage very poor - must be lots of hectads in scotland recorded recently (i.e. TP2), but not in TP1
ggsave(heavy_hecs, file="../../Maps/Upland_moth_13spp_heavy_hecs.png")

################# LOW ELEVATION BOUNDARY SHIFT - mean elevation of 10 lowest elevation hectads for each time period
## order elevation by ascending value
## take the top 10 values 
low_10_elev_heavy_hecs <- nmrsdata_heavy %>% group_by(Common_Name, time_period) %>% arrange(elevation10x10km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_heavy_hecs <- summarySE(low_10_elev_heavy_hecs, measurevar="elevation10x10km", groupvars=c("Common_Name","time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_heavy_hecs, 
          elevation10x10km[time_period == "TP1"] - elevation10x10km[time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.02098 - distribution is NOT normal
qqnorm(d)
qqline(d) ## not good

## wilcoxon signed-rank test instead
low_elev_heavy <- wilcox.test(elevation10x10km ~ time_period, data = low_elev_heavy_hecs, paired = TRUE)
low_elev_heavy ## non-significant

## plot anyway
ggpaired(low_elev_heavy_hecs, x = "time_period", y = "elevation10x10km",
         color = "time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)

## method 2: accounting for status change using regression model
# calculate average latitude of 10 lowest grid cells in each time period = range margin
# range margin change = difference in range margin between time periods (TP2 - TP1)
# positive values == shift to north
## difference in low elevation boundary between TP1 and TP2
lrb_raw_shift <- low_elev_heavy_hecs %>%
  group_by(Common_Name) %>%
  dplyr::summarise(lrb_diff = elevation10x10km[time_period == "TP2"] - elevation10x10km[time_period == "TP1"])
# negative values = srb has moved downhill
# positive values = srb has moved uphill (as predicted)

## ranges from 129m downhill (Light Knot Grass)
## to 461m uphill (Golden-rod Brindle)
## do these change after accounting for distribution change?


# distributional change = proportion of occupied hectads for each species out of all 
# hectads with observations of any species in TP1 and TP2
# change in distribution = TP2 - TP1
# no change in distribution = 0, positive = expansion, negative = contraction
all_hecs <- length(unique(nmrsdata_heavy$Hectad))
status_change1 <- nmrsdata_heavy %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(prop_hectads = n()/all_hecs)
status_change <- status_change1 %>%
  group_by(Common_Name) %>%
  dplyr::summarise(status_diff = prop_hectads[time_period == "TP2"] - prop_hectads[time_period == "TP1"])
## beech-green carpet status_diff = 0 as no change in number of occupied hectads (70) between time periods

lrb_status <- merge(lrb_raw_shift, status_change, by="Common_Name")

# regression of change in change margin ~ change in distribution
# if intercept is positive and significantly different from zero, the species group has shifting SRB towards the north
# intercept is value of SRB when status_diff is 0
hist((lrb_status$lrb_diff)^1/3) ## very right skewed - one main outlier (small autumnal moth)
lrb_status$lrb_diff_transf <- (lrb_status$lrb_diff)^1/3
## remove small autumnal and run model again - helps a bit
lrb_status <- lrb_status[!lrb_status$Common_Name=="Small Autumnal Moth",]

lrb_status_model <- lm(lrb_diff_transf ~ status_diff, data=lrb_status)
summary(lrb_status_model) 
## both non-significant

ggplot(lrb_status, aes(x=status_diff, y=lrb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in prevalence", y="Change in low elevation range margin (m)") +
  theme_classic()

# extract residuals to obtain corrected measure of SRB shifts per species 
# this is the residual shift in SRB that is not explain by the linear effect of prevalence across the studied species
# probably correlate with raw SRB estimates, but are more conservative (can show this in plot)
lrb_status_model$residuals        
# add these to srb_status dataframe
lrb_status$lrb_diff_adjust <- lrb_status_model$residuals
## plot relationship between the two
ggplot(lrb_status, aes(x=lrb_diff_adjust, y=lrb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in LRB (residuals)", y="Change in LRB (raw)") +
  theme_classic()
## very similar

## could split species into categories to show direction of responses (as per Hallfors et al 2021)
## then could do chi-squared between groups? Do species moving north also move uphill?

# srb_status$srb_category <- case_when(
#   srb_status$srb_diff_adjust>=20 ~ "Northwards",
#   srb_status$srb_diff_adjust<=-20 ~ "Southwards",
#   TRUE ~ "Stable"
# )
# proportions <- table(srb_status$srb_category)/length(srb_status$srb_category)
# percentages <- proportions*100

################# LOW LATITUDE BOUNDARY SHIFT - mean latitude of 10 most southerly hectads for each time period
## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_heavy_hecs <- nmrsdata_heavy %>% group_by(Common_Name, time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_heavy_hecs <- summarySE(low_10_lat_heavy_hecs, measurevar="lat", groupvars=c("Common_Name","time_period"))

## check for normality 
# compute the difference
d <- with(low_lat_heavy_hecs, 
          lat[time_period == "TP1"] - lat[time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # non-normal distribution
qqnorm(d)
qqline(d) ## not good

## wilcox test
low_lat_heavy <- wilcox.test(lat ~ time_period, data = low_lat_heavy_hecs, paired = TRUE)
low_lat_heavy ## significant!

## plot result
ggpaired(low_lat_heavy_hecs, x = "time_period", y = "lat",
         color = "time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
## most species moving northwards (3 southwards)

## create bar plots on elevation change over time 
ggplot(low_lat_heavy_hecs, aes(x=Common_Name, y=lat, fill=Time_period)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=lat-ci, ymax=lat+ci), width=.2, position=position_dodge(0.9))+
  xlab("Common Name") +
  ylab("Low latititude boundary") +
  theme_classic()
## not very informative - plot something else

## method 2: accounting for status change using regression model
# calculate average latitude of 10 southern-most grid cells in each time period = range margin
# range margin change = difference in range margin between time periods (TP2 - TP1)
# positive values == shift to north
low_10_northing_heavy_hecs <- nmrsdata_heavy %>% group_by(Common_Name, time_period) %>% arrange(northing) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_northing_heavy_hecs <- summarySE(low_10_northing_heavy_hecs, measurevar="northing", groupvars=c("Common_Name","time_period"))
## difference in low latitude boundary between TP1 and TP2
srb_raw_shift <- low_northing_heavy_hecs %>%
  group_by(Common_Name) %>%
  dplyr::summarise(srb_diff = northing[time_period == "TP2"] - northing[time_period == "TP1"])
# negative values = srb has moved southwards
# positive values = srb has moved northwards (as predicted)
## divide each by 1000 to give values in km - then residuals from model will also be in km
srb_raw_shift$srb_diff <- srb_raw_shift$srb_diff/1000
## ranges from 169km northwards (Pretty Pinion)
## to 105km southwards (Smoky Wave)
## do these change after accounting for distribution change?


# distributional change = proportion of occupied hectads for each species out of all 
# hectads with observations of any species in TP1 and TP2
# change in distribution = TP2 - TP1
# no change in distribution = 0, positive = expansion, negative = contraction
all_hecs <- length(unique(nmrsdata_heavy$Hectad))
status_change1 <- nmrsdata_heavy %>% group_by(Common_Name, time_period) %>%
  dplyr::summarise(prop_hectads = n()/all_hecs)
status_change <- status_change1 %>%
  group_by(Common_Name) %>%
  dplyr::summarise(status_diff = prop_hectads[time_period == "TP2"] - prop_hectads[time_period == "TP1"])

srb_status <- merge(srb_raw_shift, status_change, by="Common_Name")

# regression of change in change margin ~ change in distribution
# if intercept is positive and significantly different from zero, the species group has shifting SRB towards the north
# intercept is value of SRB when status_diff is 0
hist(srb_status$srb_diff) ## looks good
srb_status_model <- lm(srb_diff ~ status_diff, data=srb_status)
summary(srb_status_model) 
## intercept is non-significant (21km, p=0.07)
## shift in prevalence  significant (p=0.02) & negative

ggplot(srb_status, aes(x=status_diff, y=srb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in prevalence", y="Change in southern range margin (km)") +
  theme_classic()

# extract residuals to obtain corrected measure of SRB shifts per species 
# this is the residual shift in SRB that is not explain by the linear effect of prevalence across the studied species
# probably correlate with raw SRB estimates, but are more conservative (can show this in plot)
srb_status_model$residuals        
# add these to srb_status dataframe
srb_status$srb_diff_adjust <- srb_status_model$residuals
## plot relationship between the two
ggplot(srb_status, aes(x=srb_diff_adjust, y=srb_diff)) +
  geom_point() +
  geom_smooth(method = "lm", colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x="Change in SRB (residuals)", y="Change in SRB (raw)") +
  theme_classic()
## very similar

## can split species into categories to show direction of responses (as per Hallfors et al 2021)
# Northwards: shift of >= 20km
# Stable: shift beteween 20km and -20km
# Southwards: shift of <= -20km

srb_status$srb_category <- case_when(
  srb_status$srb_diff_adjust>=20 ~ "Northwards",
  srb_status$srb_diff_adjust<=-20 ~ "Southwards",
  TRUE ~ "Stable"
)
proportions <- table(srb_status$srb_category)/length(srb_status$srb_category)
percentages <- proportions*100

## put results in one dataframe 
recorded_results <- merge(mean_elev_rec, low_elev_rec, by="Common_Name", all=T)
recorded_results <- merge(recorded_results, low_lat_rec, by="Common_Name", all=T)
recorded_results$recording_level <- "Recorded"

well_rec_results <- merge(mean_elev_well, low_elev_well, by="Common_Name", all=T)
well_rec_results <- merge(well_rec_results, low_lat_well, by="Common_Name", all=T)
well_rec_results$recording_level <- "Well Recorded"

heavy_rec_results <- merge(mean_elev_heavy, low_elev_heavy, by="Common_Name", all=T)
heavy_rec_results <- merge(heavy_rec_results, low_lat_heavy, by="Common_Name", all=T)
heavy_rec_results$recording_level <- "Heavily Recorded"

nmrs_range_shift_results <- rbind(recorded_results, well_rec_results, heavy_rec_results)
write.csv(nmrs_range_shift_results, file="../../Results/nmrs_range_shift_results.csv", row.names=FALSE)









## now calculate each species prevalence for T1 and T2
## then subtract the SRB (southern range boundary, i.e. average latitude of 10 southernmost hectads) in T2 from T1 = change in SRB
## then the change in prevalence by subtracting prevalence in T2 from T1
## then model change in km between TPs as a function of change in prevalence
## change in SRB ~ change in prevalence 
## if intercept is positive and significantly different from zero, the species group has shifted their SRB


## for now calculate centroid of each species distribution - average lat and lon
## for each time period, then calculate the bearings between the two
## this gives us an idea of what direction a species range is moving

## centroid shift
### TP 1
species <- nmrsdata_tp1_rec$Common_Name
long <- nmrsdata_tp1_rec$lon
lat <- nmrsdata_tp1_rec$lat
data <- as.data.frame(cbind(long,lat,species))

data$lat <- as.numeric(data$lat)
data$long <- as.numeric(data$long)

library(geosphere)
library(ggplot2)
library(dplyr)
cntrd <- function(x) {
  data.frame(centroid(as.matrix(x[,c("long", "lat")])))
}

by(data, data$species, cntrd) %>% head()
centroids_geo_tp1 <- group_by(data, species) %>%
  do(cntrd(.))

## TP2
species <- nmrsdata_tp2_rec$Common_Name
long <- nmrsdata_tp2_rec$lon
lat <- nmrsdata_tp2_rec$lat
data <- as.data.frame(cbind(long,lat,species))

data$lat <- as.numeric(data$lat)
data$long <- as.numeric(data$long)

cntrd <- function(x) {
  data.frame(centroid(as.matrix(x[,c("long", "lat")])))
}

by(data, data$species, cntrd) %>% head()
centroids_geo_tp2 <- group_by(data, species) %>%
  do(cntrd(.))

tp1_cent <- centroids_geo_tp1[,c(2,3)]
tp2_cent <- centroids_geo_tp2[,c(2,3)]

cent_rec <- merge(centroids_geo_tp1, centroids_geo_tp2, by="species")
x <- bearing(tp1_cent, tp2_cent)
cent_rec$bearing <- x

dir <- setNames( seq(0, 315 , by=45),
                 c("North", "North-east", "East" ,
                   "South-east", "South", "South-west", 
                   "West", "North-west"))

vec = sample(  c("North", "North-east", "East" ,
                 "South-east", "South", "South-west", 
                 "West", "North-west"), 25, TRUE)

dir[vec]


nms <-  c("North", "North-east", "East" ,
          "South-east", "South", "South-west", 
          "West", "North-west")
dir <- seq(0, 315 , by=45)
mt <- match(vec, nms) 

course <- (x + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
cent_rec$BEARING <- x
cent_rec$COURSE  <- course
cent_rec

d2c.2 <- function(x) {
  upper <- seq(from = 11.25, by = 22.5, length.out = 17)
  card1 <- c('N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', 'N')
  ifelse(x>360 | x<0,NA,card1[findInterval(x,upper,rightmost.closed = T)+1])
}

cent_rec$direction <- d2c.2(cent_rec$COURSE)

##




library(spocc)
# Get occurrences for Phyllomedusa
occurrrences <- occ(query = 'Phyllomedusa', from = 'gbif', limit = 5000)
data <- occurrrences$gbif$data$Phyllomedusa
# Remove NAs
remove_na <- is.na(data$longitude) | is.na(data$latitude)
data <- data[!remove_na, ]
xy <- data[, c("longitude", "latitude")] # coordinates
sp_name <- data$name







library(sf)

library(rgdal)
gb_10km <- readOGR(dsn="../../Data/gb_10km.shp")
## xmin: 2780000 ymin: 2880000 xmax: 3930000 ymax: 4590000
gb_10km@proj4string
gb_10km <- spTransform(gb_10km, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


st_crs(gb_10km)
ggplot() + 
  geom_sf(data = gb_10km, size = 3, color = "black", fill = "cyan1") + 
  ggtitle("AOI Boundary Plot") + 
  coord_sf()


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

ConvertCoordinates(2780000, 2880000)

####################################################################################################
rm(list = ls())
## subset species by habitat association - species coded as 1 or 2 for upland habitat (obligate use and predominate use respectively)

## read in data
trait_data <- read.csv("../../Data/ecological_traits.csv", header=TRUE)
nmrsdata <- readRDS("../../Data/NMRS/NMRS_1km_elevation_temp.rds")

length(unique(nmrsdata$Common_Name)) ## 822 species 
length(unique(trait_data$common_name)) ## 968 species 

## subset species by upland levels 1 and 2 
north_species <- subset(trait_data, X7_montane_upland == 1 | X7_montane_upland == 2)
## remove nymphalids (butterflies)
north_species <- north_species[!(north_species$family=="Nymphalidae"),] ## 36 species
north_species <- north_species[,c(3,4,179)]
#write.csv(north_species, "../../Data/north_species2.csv", row.names=F)

## subset species in nmrs data
nmrsdata <- merge(nmrsdata, north_species, by.x="Preferred_Name", by.y="scientific_name")
length(unique(nmrsdata$Common_Name)) ## 35 species
## Small Lappet not in NMRS data - only records are pre-1970

## calc mean elevation and temp for each species over time & space
nmrs_mean_elev_temp <- nmrsdata %>% group_by(Common_Name) %>% summarize(mean_elev = mean(elevation1x1km),
                                                                        mean_temp = mean(temperature))
range(nmrs_mean_elev_temp$mean_elev) ## 350m - 1289m
range(nmrs_mean_elev_temp$mean_temp) ## 4.5C - 9.56C

## yellow-ringed carpet and manchester treble-bar (sub-arctic moths not included in above list)
## they fall within both ranges above but are classed as upland==3

worldmap = map_data('world')
for(i in unique(nmrsdata$Preferred_Name)) {
  print(i)
  # Printing ggplot within for-loop
  
  temp_plot <- ggplot() + 
    geom_polygon(data = worldmap, 
                 aes(x = long, y = lat, group = group), 
                 fill = 'gray90', color = 'black') + 
    coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
    theme_void() + 
    geom_point(data = nmrsdata[nmrsdata$Preferred_Name==i,], 
               aes(x = as.numeric(lon), 
                   y = as.numeric(lat), colour=elevation10x10km), size=1) +
    scale_color_viridis_c(name="Elevation") + 
    ggtitle(i) +
    theme(title = element_text(size = 12))
  
  ggsave(temp_plot, file=paste0("../../Maps/NMRS_10km_elevation_", i,".png"), width = 20, height = 25, units = "cm")
  Sys.sleep(2)
}

## save nmrs data
write.csv(nmrsdata, file="../../Data/nmrsdata_northern_spp_10km.csv", row.names=FALSE)
## 4,659,318 rows

## calc mean elevation for each species - get an idea of if we will detect elevational trends
library(dplyr)
nmrs_mean_elev <- nmrsdata %>% group_by(Preferred_Name) %>% summarize(mean_elev = mean(elevation1x1_new))
north_species <- merge(nmrs_mean_elev, north_species, by.x="Preferred_Name", by.y="Scientific_Name")
north_species <- north_species[,c(3,1,2)]
write.csv(north_species, "../../Data/north_elevation_1km.csv", row.names=FALSE)
## minimum elevation across years and 10km grids = 201m (spp = Epirrhoe galiata)
## therefore we should capture elevation shifts well


##################################################################    
rm(list = ls())

options(scipen=999)

library(raster)
library(ncdf4)
library(rgdal) # package for geospatial analysis
##################################################################             
## open gridded climate data from HadUK 1x1km 30        
ncpath <- "../../Data/"
ncname <- "tas_hadukgrid_uk_1km_ann-30y_198101-201012"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "tas"  

tmp_raster <- brick(ncfname, varname="tas")
tmp_raster; class(tmp_raster)

tmp_raster[tmp_raster %in% 9.969209968386869047442886268468442020e+36] <- NA 

plot(tmp_raster)
res(tmp_raster)
crs(tmp_raster)
r_newproj <- projectRaster(tmp_raster,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(r_newproj, axes = T)
average_temp <- as.data.frame(rasterToPoints(r_newproj)) ## get temp values as dataframe
colnames(average_temp) <- c("lon", "lat", "temperature")

## better plot of temperature data
library(ggplot2)
worldmap = map_data('world')
temp_map <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = average_temp, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour = temperature), size=1) + 
  scale_color_viridis_c(name="Mean air temperature") + 
  theme_void() +
  theme(title = element_text(size = 12))
temp_map

## plot temperature of all NMRS 1km sites
nmrsdata_1km <- readRDS("../../Data/NMRS/NMRS_1km_cleaned.rds")
lat_lon <- nmrsdata_1km[,c(8:9)]
lat_lon <- unique(lat_lon)
rm(nmrsdata_1km)
rm(r_newproj)
rm(tmp_raster)

temp_latlon <- average_temp[,1:2]

## need to source Find Closest Rcpp function first!

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


result <- find_closest(lat_lon[, 2], lat_lon[, 1], 
                       temp_latlon[, 2], temp_latlon[, 1])

result ### works super fast!!! assume numbers are the corresponding row numbers that each lat_lon value matches in temp_latlon?
## 62539 results, but 59967 unique values - some temperature sites have matched with multiple NMRS sites (mainly on the coast & Spurn Head)
## add row numbers to temp_latlon
lat_lon$row_numbers <- rownames(lat_lon)
lat_lon$temp_row <- result
average_temp$row_numbers <- rownames(average_temp)

lat_lon_final <- merge(lat_lon, average_temp, by.x="temp_row", by.y="row_numbers")
lat_lon_final <- lat_lon_final[,-c(1,4:6)]
colnames(lat_lon_final) <- c("lon", "lat", "temperature")


## map of mean air temp at NMRS 1km sqaures/points
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = lat_lon_final, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour = temperature), size=1) + 
  scale_color_viridis_c(name="30 year average temperature") + 
  theme_void() +
  theme(title = element_text(size = 12))

## save temp data
write.csv(lat_lon_final, file="../../Data/1km_NMRS_30y_mean_temp.csv", row.names=FALSE)

#### merge temperature data back into nmrs data
nmrsdata_1km <- readRDS("../../Data/NMRS/nmrsdata_elevation_1km.rds")
temp_data <- read.csv("../../Data/1km_NMRS_30y_mean_temp.csv", header=TRUE)

nmrsdata_1km <- merge(nmrsdata_1km, temp_data, by=c("lat", "lon"))
## save file
saveRDS(nmrsdata_1km, "../../Data/NMRS/NMRS_1km_elevation_temp.rds")


###################################################################################################
nmrsdata_1km <- read.csv("../../Data/NMRS/NMRS_1km_elevation.csv", header=TRUE)

## calculate mean and SD of temperature in each species range
## using for loop

## take out year info from nmrs data? temp and elevation does not change over time so not relevant
library(plyr)
colnames(nmrsdata_1km)[3] <- "Grid_Square"
nmrsdata_species <- ddply(nmrsdata_1km, .(Grid_Square, Common_Name, lat, lon), numcolwise(mean)) ## removes species data, just hectads recorded for each year
## takes a while

climate_envelope <- NULL
species <- unique(nmrsdata_species$Common_Name) ## 875 species

start_time <- Sys.time()
for (i in species) {
  species_data <- nmrsdata_species[nmrsdata_species$Common_Name==i,] ## species of interest
  mean_temp <- mean(species_data$temperature)
  min_temp <- min(species_data$temperature)
  max_temp <- max(species_data$temperature)
  sd_temp <- sd(species_data$temperature)
  mean_elev <- mean(species_data$elevation1x1km)
  climate_tmp <- data.frame(i, mean_temp, sd_temp, min_temp, max_temp, mean_elev)
  climate_envelope <- rbind(climate_tmp, climate_envelope)
}## loop through each species
## subset data by species
## calc mean, sd and range of temperature
## collate data
end_time <- Sys.time()
end_time - start_time ## ~30 secs

## calc 33% of max. mean temperature
## take species up that 20% value
## e.g. if max = 12, 20% = 2.4
max(climate_envelope$mean_temp)*0.333
## 4.042917
min(climate_envelope$mean_temp) + 4.042917
## 8.298319

## include species in the lower quartile of data (25%)
upland_climate_envelope <- climate_envelope[climate_envelope$mean_temp <= quantile(climate_envelope$mean_temp , 0.2 ),]
upland_climate_envelope <- climate_envelope2[climate_envelope2$mean_elev>=300,] 
## 245 species (can limit it to 200, but after temperature filter of 8.5C 
## species between 200 and 300m are removed anyway)
upland_climate_envelope <- upland_climate_envelope[upland_climate_envelope$mean_temp<=8.29,] ## 44 species


library(dplyr)
climate_envelope <- climate_envelope %>% arrange(mean_temp)
obs <- nrow(climate_envelope) 
climate_envelope2 <- climate_envelope %>% filter(row_number() < obs * 0.1)


quantile(climate_envelope$mean_temp, 0.05)
