##########################
#### user: Lisbeth Hordley
#### date: August 2021
#### info: Range shift analysis/exploration using BNM data

rm(list = ls())

## load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(sf)
library(ggplot2)
library(plyr)
library(rnrfa)

## read in BNM data
bnmdata <- readRDS("../../Data/BNM/BNM_UK_data_1970-2019_processed_minimalcol.rds")

## remove Northern Ireland data by Country
bnmdata<-bnmdata[!(bnmdata$Country=="Northern Ireland"),]
## remove unecessary columns 
bnmdata[,c("Date","Week","Month","Country","listL"):=NULL] 
range(bnmdata$Year) ## 1970 - 2019 
## different time intervals? or stick to same as moths? 

## get lat/lon and easting/northing at 10km resolution using 10k grid ref
colnames(bnmdata)[2] <- "Hectad"

## convert gridref to east/north - gives easting and northing at 10km scale - used to measure distances between sites later
gridrefs <- unique(as.vector(bnmdata[,2])) ## 2786 10km gridrefs
lon_lat <- as.data.frame(osg_parse(grid_refs=gridrefs$Hectad, coord_system = "WGS84"))
east_north <- as.data.frame(osg_parse(grid_refs=gridrefs$Hectad))
gridref_latlon_eastnorth <- cbind(gridrefs, lon_lat, east_north)
## merge back into bnm data
bnmdata <- merge(bnmdata, gridref_latlon_eastnorth, by="Hectad")
head(bnmdata)
## remove other easting and northing (which are at 1km scale)
bnmdata[,c("EASTING","NORTHING","SQ1km"):=NULL] 

## is there any missing data?
bnmdata <- bnmdata %>%
  mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))
miss_data <- bnmdata[rowSums(is.na(bnmdata)) > 0, ] 
## nothing

## add in elevation data at 10km scale
library(raster)
tmp = raster("../../Data/elevation1x1_new.tif")
crs(tmp) ## projection
res(tmp) ## resolution (1000 x 1000m) (needs multiplying by a factor of 10 to be 10,000 x 10,000m)
tmp_agg <- aggregate(tmp, fact=10)
res(tmp_agg) ## 10000 x 10000m = 10x10km
plot(tmp_agg)
hist(values(tmp_agg))

## re-project
tmp_newproj <- projectRaster(tmp_agg,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(tmp_newproj, axes = T)
tmp_elev <- raster::extract(tmp_newproj, gridref_latlon_eastnorth[,2:3], df=T)
tmp_elev <- cbind(tmp_elev, gridref_latlon_eastnorth)
## multiply elevation by 10 - see metadata https://www.eea.europa.eu/data-and-maps/data/digital-elevation-model-of-europe?fbclid=IwAR1yqtoCIjiV1QIe30h8DCaYhWH-7GO0X68CKLAkD5FY9e5WMs6OPU9YSas 
tmp_elev$elevation1x1_new <- tmp_elev$elevation1x1_new*10 ## elevation in meters 
tmp_elev <- tmp_elev[,-1]
colnames(tmp_elev)[1] <- "elevation10x10km"

## add into BNM data with 4 species
## only need hectad info to merge
tmp_elev <- tmp_elev[,-c(3:6)]
## if throwing memory issues: memory.limit(size = 15000)
bnmdata <- merge(bnmdata, tmp_elev, by="Hectad", all=T)
## save elevation data
saveRDS(bnmdata, file="../../Data/BNM/BNMdata_10km_elevation.rds")

rm(list = ls())

bnmdata <- readRDS("../../Data/BNM/BNMdata_10km_elevation.rds")

### explore data to determine which time period to use 
# count the number of recorded hectads and the number of species*hectad records in each year 
library(plyr)
bnmdata$PRESENCE <- 1
#colnames(nmrsdata)[1] <- "Grid_Square"
bnm_rec_hectad_years <- ddply(bnmdata, .(Hectad, Year, easting, northing), numcolwise(mean)) ## removes species data, just hectads recorded for each year

recording <- data.frame(Year = numeric(),
                        Hectad = numeric(),
                        Records = numeric())

for (n in years){
  year_hec <- bnm_rec_hectad_years[which(bnm_rec_hectad_years$Year==n), ] ## NMRS data for year of interest (n)
  Hectad <- nrow(year_hec) ## number of grid squares recorded in year of interest
  
  year_rec <- bnmdata[which(bnmdata$Year==n), ] ## find records for unique year sites in nmrs data (i.e. with species data)
  Records <- nrow(year_rec) ## number of records for year of interest
  
  Year <- n
  
  out <- data.frame(cbind(Year,Hectad,Records))
  recording <- rbind(recording,out)
  
  
}

summary(recording)


# now plot the number of hectads and records over time
plot(Hectad ~ Year, data = recording)
plot(Records ~ Year, data = recording)
total_hecs <- max(recording$Hectad) # 2445 (in 2014)

recording$Perc <- recording$Hectad*100/total_hecs

plot(Perc ~ Year, data = recording)


## need to decide on time periods to investigate
## current vs pre-2000
## following hickling 2006 - 25 year period between mid-points of time periods
## e.g. current is 2012-2016 (mid-point is 2014)
## 1989 would be previous mid-point
## other time period 1987-1991
## longer first time period? By the time we look at heavily recorded hectads, we might have barely any
## hectads left in first TP

## follow same as moths for now

tp1 <- 1975:1991
tp2 <- 2012:2016

# first follow Hickling et al directly by calculating which grid squares were recorded at least once in each interval
# first select out the intervals
bnm_rec_hectad_years_tp1 <- bnm_rec_hectad_years[which(bnm_rec_hectad_years$Year %in% tp1), ]
bnm_rec_hectad_years_tp2 <- bnm_rec_hectad_years[which(bnm_rec_hectad_years$Year %in% tp2), ]

# then get rid of the year information to leave a non-redundant list of grid squares recorded in each interval
bnm_rec_hectad_tp1 <- ddply(bnm_rec_hectad_years_tp1, .(Hectad), numcolwise(mean))
bnm_rec_hectad_tp2 <- ddply(bnm_rec_hectad_years_tp2, .(Hectad), numcolwise(mean))
## this just takes the mean year for each grid square - just gives a unique list of hectads recorded in each time period
## 2127 in tp1
## 2702 in tp2

# plot these on a map

ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =bnm_rec_hectad_tp1,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat)), size=1, colour="blue") +
  geom_point(data =bnm_rec_hectad_tp2,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat)), size=1, colour="red") +
  theme(title = element_text(size = 12))
## good coverage - very few hectads not covered

## now we want to recombine these datasets in a way that ditches grid squares not present in both TPs
bnm_rec_hectad_tp1_vec <- as.character(bnm_rec_hectad_tp1$Hectad)
bnm_rec_hectad_tp2_vec <- as.character(bnm_rec_hectad_tp2$Hectad)

bnm_rec_hectad_full <- append(bnm_rec_hectad_tp1_vec,bnm_rec_hectad_tp2_vec) 

bnm_rec_hectad_full <- data.frame(cbind(bnm_rec_hectad_full,1)) ## list of all hectads in both TPs
colnames(bnm_rec_hectad_full) <- c("Hectad","Intervals")
bnm_rec_hectad_full$Hectad <- as.factor(as.character(bnm_rec_hectad_full$Hectad)) ## 4829 (includes duplicates)
bnm_rec_hectad_full$Intervals <- as.numeric(as.character(bnm_rec_hectad_full$Intervals))


bnm_rec_hectad_good <- ddply(bnm_rec_hectad_full, .(Hectad), numcolwise(sum)) ## 2723 unique hectads ## hectads either have 1 (only recorded in one TP), or 2 (recorded in both TPs)

bnm_rec_hectad_good <- bnm_rec_hectad_good[which(bnm_rec_hectad_good$Intervals > 1), ] ## only keep grid squares recorded in both (i.e. intervals > 1)
## 2106 good hectads

# drop unused levels
bnm_rec_hectad_good <- droplevels(bnm_rec_hectad_good)

# generate the vector
bnm_good_hectads <- levels(bnm_rec_hectad_good$Hectad)

# pick out the good records
length(unique(bnmdata$Hectad)) ## 2786 total hectads
bnmdata_good_hectads <- bnmdata[which(bnmdata$Hectad %in% bnm_good_hectads), ]
length(unique(bnmdata_good_hectads$Hectad)) ## 2106
# we lose 24% of hectads

bnm_good_hectads2 <- ddply(bnmdata_good_hectads, .(Hectad, Year), numcolwise(mean))

## plot all 2106 hectads (all hectads which are recorded in both time periods)
plot <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =bnm_good_hectads2, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), size=1, colour="purple") +
  theme(title = element_text(size = 12))
plot


## we can also further restrict to the hectads that have records of a minimum percentage of the regional species richness
# this is going to take quite a complex loop to construct:
# for every hectad in the NMRS_GB_good dataframe, we need to calculate the observed species richness
# then we need to calculate pairwise distances to all other recorded hectads, pick out the nearest 100, 
# and calculate the regional species richness of these 100
# and finally turn this into a percentage of regional SR that is recorded in the hectad

# start with the list of recorded hectads, re-attach the eastings and northings
bnm_good <- ddply(bnm_rec_hectad_years, .(Hectad,easting,northing), summarise,
                   COUNT = sum(PRESENCE))

# select out only the good hectads
bnm_good <- bnm_good[which(bnm_good$Hectad %in% bnm_good_hectads), ] ## 2106 hectads
time_period1 <- as.data.frame(tp1)
colnames(time_period1)[1] <- "year"
time_period1$interval <- "tp1"
time_period2 <- as.data.frame(tp2)
colnames(time_period2)[1] <- "year"
time_period2$interval <- "tp2"
time_period <- rbind(time_period1, time_period2)

bnm_good$easting <- as.numeric(bnm_good$easting)
bnm_good$northing <- as.numeric(bnm_good$northing)

# now start to construct the loop:
hec_records <- data.frame(Year = numeric(),
                          Hectad = factor(),
                          Easting = numeric(),
                          Northing = numeric(),
                          hec_SR = numeric(),
                          reg_SR = numeric(),
                          hec_perc = numeric())
i <- 1

for (x in bnm_good_hectads){
  if (round(i, -1) == i){
    print(i)
  }
  hec <- bnm_good[which(bnm_good$Hectad == x), ]    # pick out details of focal hectad
  candidates <- bnm_good[which(bnm_good$Hectad != x), ]   # pick out details of all others
  
  hec_east <- hec[[1,2]]     # easting of target
  hec_north <- hec[[1,3]]    # northing of target
  
  candidates$east_diff <- candidates$easting - hec_east             # longitudinal difference
  candidates$north_diff <- candidates$northing - hec_north          # latitudinal difference
  
  candidates$distance <- sqrt((candidates$east_diff^2) + (candidates$north_diff^2))    # absolute difference
  
  candidates <- candidates[order(candidates$distance),] # sort by distance ascending
  
  closest <- candidates[1:100,] # select out closest 100
  
  # now we want to calculate species richness from the hectad and from the region in both time periods
  for (n in unique(time_period$interval)){
    y <- time_period$year[time_period$interval==n] ## tp years of interest
    year_recs <- bnmdata_good_hectads[bnmdata_good_hectads$Year %in% y, ] ## nmrs good hectad data for tp of interest
    
    hec_recs <- year_recs[which(year_recs$Hectad == x), ] # pull out hectad records
    hec_SR <- length(unique(hec_recs$Species))       # calculate species richness
    
    reg_recs <- year_recs[which(year_recs$Hectad %in% closest$Hectad), ]  # pull out region records
    reg_recs <- rbind(reg_recs,hec_recs)    # add in hectad records (they're part of the regional richness too!)
    reg_SR <- length(unique(reg_recs$Species))
    
    hec_perc <- hec_SR*100/reg_SR
    
    out <- cbind(x,n,hec_east,hec_north,hec_SR,reg_SR,hec_perc)
    hec_records <- rbind(hec_records,out)
  }    
  i <- i+1
}

colnames(hec_records) <- c("HECTAD","TIME.PERIOD","EASTING","NORTHING","HECTAD.SR","REGION.SR","PERCENT.RECORDED")
summary(hec_records)

# make things that should be numeric, numeric
hec_records$EASTING <- as.numeric(as.character(hec_records$EASTING))
hec_records$NORTHING <- as.numeric(as.character(hec_records$NORTHING))
hec_records$HECTAD.SR <- as.numeric(as.character(hec_records$HECTAD.SR))
hec_records$REGION.SR <- as.numeric(as.character(hec_records$REGION.SR))
hec_records$PERCENT.RECORDED <- as.numeric(as.character(hec_records$PERCENT.RECORDED))

summary(hec_records)

# most of the hectad*year combinations are "well recorded" but only slightly under a quarter are "heavily recorded"

# label each hectad according to its recording level in each time period
hec_records$RECORDING.LEVEL <- as.factor(ifelse(hec_records$PERCENT.RECORDED >= 25, "Heavily recorded",
                                                ifelse(hec_records$PERCENT.RECORDED >= 10, "Well recorded",
                                                       ifelse(hec_records$PERCENT.RECORDED >0, "Recorded",
                                                              "Not recorded"))))
## get lat/lon values again to plot with ggplot
lat_lon <- bnmdata_good_hectads %>% distinct(Hectad, lat, lon, .keep_all = FALSE)

hec_records <- merge(hec_records, lat_lon, by.x="HECTAD", by.y="Hectad")

# this data takes a while to generate so let's back it up
write.csv(hec_records, "../../Data/BNM/Hectad_recording_levels_1975_1991_2012_2016.csv", row.names = F)

hec_records <- read.csv("../../Data/BNM/Hectad_recording_levels_1975_1991_2012_2016.csv", header=TRUE)

## Analysis 1: hectads with at least one species recorded in BOTH time periods: 2106 hectads
## file = bnmdata_good_hectads (needs subsetting by years)
years <- c(1975:1991, 2012:2016)
bnmdata_rec_hecs <- bnmdata_good_hectads[which(bnmdata_good_hectads$Year %in% years), ]
length(unique(bnmdata_rec_hecs$Hectad)) ## 2106 hectads
length(unique(bnmdata_rec_hecs$Year)) ## 22 years
## save file
saveRDS(bnmdata_rec_hecs, "../../Data/BNM/BNMdata_rec_hecs_75_91_12_16.rds")

## Analysis 2a: hectads classed as well-recorded in BOTH time periods: 41 hectads (not doable)
well_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Well recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
well_hecs <- well_hecs[which(well_hecs$n_row > 1), ] 
well_hecs$n_row <- NULL

## Analysis 2b: hectads classed as heavily recorded in BOTH time periods: 1209 hectads
heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
heavy_hecs <- heavy_hecs[which(heavy_hecs$n_row > 1), ] 
heavy_hecs$n_row <- NULL
## merge into nmrs data and filter by years so only have years within time periods
bnmdata_heavy_hecs <- bnmdata[which(bnmdata$Hectad %in% heavy_hecs$HECTAD)]
length(unique(bnmdata_heavy_hecs$Hectad)) ## 1209 hectads
bnmdata_heavy_hecs <- bnmdata_heavy_hecs[which(bnmdata_heavy_hecs$Year %in% years)]
length(unique(bnmdata_heavy_hecs$Year)) ## 22 years
## save file
saveRDS(bnmdata_heavy_hecs, "../../Data/BNM/BNMdata_heavy_hecs_75_91_12_16.rds")

## Analysis 2c: hectads classed as well OR heavily recorded in BOTH time periods: 1651 hectads
well_heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded" | RECORDING.LEVEL=="Well recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
well_heavy_hecs <- well_heavy_hecs[which(well_heavy_hecs$n_row > 1), ] 
well_heavy_hecs$n_row <- NULL
## merge into nmrs data and filter by years so only have years within time periods
bnmdata_well_heavy_hecs <- bnmdata[which(bnmdata$Hectad %in% well_heavy_hecs$HECTAD)]
length(unique(bnmdata_well_heavy_hecs$Hectad)) ## 1651 hectads
bnmdata_well_heavy_hecs <- bnmdata_well_heavy_hecs[which(bnmdata_well_heavy_hecs$Year %in% years)]
length(unique(bnmdata_well_heavy_hecs$Year)) ## 22 years
## save file
saveRDS(bnmdata_well_heavy_hecs, "../../Data/BNM/BNMdata_well_heavy_hecs_75_91_12_16.rds")

###############################################################
## recorded hectads
rm(list = ls())

bnmdata_rec_hecs <- readRDS("../../Data/BNM/BNMdata_rec_hecs_75_91_12_16.rds")
## 2106 hectads in TP1 and TP2

## Filter 1: northern butterflies only: mountain ringlet, scotch argus, northern brown argus and large heath
northern_butterflies <- c("Large Heath", "Scotch Argus", "Mountain Ringlet", "Northern Brown Argus")
bnmdata_rec_hecs <- bnmdata_rec_hecs[which(bnmdata_rec_hecs$Species %in% northern_butterflies), ]
length(unique(bnmdata_rec_hecs$Species)) ## 4 species
length(unique(bnmdata_rec_hecs$Hectad)) ## 472 hectads

####### Filter 2: species must have enough coverage to calculate range shifts
####### species need to occupy at least 20 hectads in each time period
bnmdata_rec_hecs$Time_period <- ifelse(bnmdata_rec_hecs$Year==2012 | bnmdata_rec_hecs$Year==2013 | bnmdata_rec_hecs$Year==2014
                                        | bnmdata_rec_hecs$Year==2015 | bnmdata_rec_hecs$Year==2016, "TP2", "TP1")
# get rid of the year information, keeping only which hectads each species was recorded in:
bnmdata_rec_hecs2 <- ddply(bnmdata_rec_hecs, .(Hectad,Species,easting,northing,Time_period), summarise,
                            PRESENCE = 1) 
# and now count the number of hectads per species for each year
tp_species_hecs <- ddply(bnmdata_rec_hecs2, .(Species, Time_period), summarise,
                         HECTADS = sum(PRESENCE))
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_hecs <- tp_species_hecs %>% group_by(Species) %>% filter(all(HECTADS>=20))
species_keep <- unique(tp_species_hecs$Species) ## all 4 species still
## remove duplicated records at 10km resolution - same species/year/hectad (probably have different 1km sites)
bnmdata_rec_hecs <- unique(bnmdata_rec_hecs) ## 1842

## map
worldmap = map_data('world')
rec_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =bnmdata_rec_hecs, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
rec_hecs
## save graph
ggsave(rec_hecs, file="../../Maps/Northern_butterflies_4spp_rec_hecs.png")

### now calculate centroid elevation shift between time periods
## centroid = average elevation of all sampled occurrences during each time period
## and measure the low elevation shift - mean elevation of 10 lowest elevation hectads for each time period

## mean elevation for each species and each time period + SE + confidence intervals
library(Rmisc)
mean_elev_rec_hecs <- summarySE(bnmdata_rec_hecs, measurevar="elevation10x10km", groupvars=c("Species","Time_period"))

## check for normality 
# compute the difference
d <- with(mean_elev_rec_hecs, 
          elevation10x10km[Time_period == "TP1"] - elevation10x10km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.8848 - distribution of the differences are not significantly different from normal distribution 
## i.e. they follow a normal distribution so t-test can be used
qqnorm(d)
qqline(d) ## looks fine

mean_elev_rec <- t.test(elevation10x10km ~ Time_period, data = mean_elev_rec_hecs, paired = TRUE)
mean_elev_rec ## non-significant

## plot anyway
ggpaired(mean_elev_rec_hecs, x = "Time_period", y = "elevation10x10km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Species")+
  xlab("Time period")+
  ylab("Mean elevation")+
  stat_compare_means(method="t.test", paired = TRUE)
## 3 out of 4 species moving downhill
## probably not enough statistical power to capture significant change overall

## create bar plots on elevation change over time 
ggplot(mean_elev_rec_hecs, aes(x=Species, y=elevation10x10km, fill=Time_period)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=elevation10x10km-ci, ymax=elevation10x10km+ci), width=.2, position=position_dodge(0.9))+
  xlab("Common Name") +
  ylab("Mean elevation") +
  theme_classic()

#####################################
## low elevational boundary shift
## this is the 10 lowest OCCUPIED hectads - can't have the same hectad repeated in this 10 just because it was recorded in multiple years

## first remove year info - this means we can only look at occupied hectads within each time period, rather than each year
## e.g. if one hectad is recorded in multiple years within one time period, it would be included multiple times in the 10 lowest hectads
bnmdata_rec_hecs2 <- ddply(bnmdata_rec_hecs, .(Hectad,Species,easting,northing,Time_period,elevation10x10km,lat,lon), summarise,
                           PRESENCE = 1) 
## now order elevation by ascending value
## take the top 10 values 
low_10_elev_rec_hecs <- bnmdata_rec_hecs2 %>% group_by(Species, Time_period) %>% arrange(elevation10x10km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_rec_hecs <- summarySE(low_10_elev_rec_hecs, measurevar="elevation10x10km", groupvars=c("Species","Time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_rec_hecs, 
          elevation10x10km[Time_period == "TP1"] - elevation10x10km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.5566 - distribution is normal
qqnorm(d)
qqline(d) ## looks fine

## paired t-test
low_elev_rec <- t.test(elevation10x10km ~ Time_period, data = low_elev_rec_hecs, paired = TRUE)
low_elev_rec ## non-significant (0.07)

## plot anyway
ggpaired(low_elev_rec_hecs, x = "Time_period", y = "elevation10x10km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Species")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
## all 4 species expand their distribution downhill

## create bar plots on elevation change over time 
ggplot(low_elev_rec_hecs, aes(x=Species, y=elevation10x10km, fill=Time_period)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=elevation10x10km-ci, ymax=elevation10x10km+ci), width=.2, position=position_dodge(0.9))+
  xlab("Common Name") +
  ylab("Low elevation boundary") +
  theme_classic()

########################################
## Latitudinal shifts (trailing edge boundary)
## find 10 most southerly hectads & calculate mean latitude across years

## use same data frame as above that has year info removed
## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_rec_hecs <- bnmdata_rec_hecs2 %>% group_by(Species, Time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_rec_hecs <- summarySE(low_10_lat_rec_hecs, measurevar="lat", groupvars=c("Species","Time_period"))
## higher latitude = more northerly
## mountain ringlet does move northerly, but not to higher elevation... (maybe it does, just can't see at 10km scale)

## check for normality 
# compute the difference
d <- with(low_lat_rec_hecs, 
          lat[Time_period == "TP1"] - lat[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.3528 - normal distribution
qqnorm(d)
qqline(d) ## looks good

## paired t-test
low_lat_rec <- t.test(lat ~ Time_period, data = low_lat_rec_hecs, paired = TRUE)
low_lat_rec ## non-significant

## plot anyway
ggpaired(low_lat_rec_hecs, x = "Time_period", y = "lat",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Species")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="t.test", paired = TRUE)
## 2 move north, 2 move south

## create bar plots on elevation change over time 
ggplot(low_lat_rec_hecs, aes(x=Species, y=lat, fill=Time_period)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=lat-ci, ymax=lat+ci), width=.2, position=position_dodge(0.9))+
  xlab("Common Name") +
  ylab("Low latititude boundary") +
  theme_classic()

#############################################
## direction of centroid shift 

bnmdata_rec_hecs_tp1 <- bnmdata_rec_hecs[bnmdata_rec_hecs$Time_period=="TP1",]
bnmdata_rec_hecs_tp2 <- bnmdata_rec_hecs[bnmdata_rec_hecs$Time_period=="TP2",]

## first time period
species <- bnmdata_rec_hecs_tp1$Species
long <- bnmdata_rec_hecs_tp1$lon
lat <- bnmdata_rec_hecs_tp1$lat
data <- as.data.frame(cbind(long,lat,species))

data$lat <- as.numeric(data$lat)
data$long <- as.numeric(data$long)

library(letsR)
PAM <- lets.presab.points(data[,1:2], data[,3], xmn=min(data$long), xmx=max(data$long),
                          ymn=min(data$lat), ymx=max(data$lat))
plot(PAM)
plot(PAM, name="Mountain Ringlet")
summary(PAM)
centroids_tp1_rec <- lets.midpoint(PAM, planar=FALSE, method="MCC")
colnames(centroids_tp1_rec) <- c("species","long","lat")
centroids_tp1_rec$cent <- "yes"
data$cent <- "no"

final_centroid <- rbind(data,centroids_tp1_rec)
species <- unique(final_centroid$species)
## plot map for each species with red dot signifying centroid
## colour=cent 
## save graphs to look at
worldmap = map_data('world')

for(i in species) {
  print(i)
  # Printing ggplot within for-loop
  
  temp_plot <- ggplot() + 
    geom_polygon(data = worldmap, 
                 aes(x = long, y = lat, group = group), 
                 fill = 'gray90', color = 'black') + 
    coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
    theme_void() + 
    geom_point(data = final_centroid[final_centroid$species==i,], 
               aes(x = as.numeric(long), 
                   y = as.numeric(lat), colour=cent), size=2) +
    scale_color_manual(values=c("black", "red")) +
    labs(colour = "Centroid value") 
  ggtitle(i) +
    theme(title = element_text(size = 12))
  
  ggsave(temp_plot, file=paste0("../../Maps/Centroids/BNM_centroid_TP1_", i,".png"), width = 20, height = 25, units = "cm")
  Sys.sleep(2)
}
## all look reasonable - go with this for now
## once we've got direction for each species, plot centroid for each TP and check direction
## use geom_line to connect the two dots

## second time period
species <- bnmdata_rec_hecs_tp2$Species
long <- bnmdata_rec_hecs_tp2$lon
lat <- bnmdata_rec_hecs_tp2$lat
data <- as.data.frame(cbind(long,lat,species))

data$lat <- as.numeric(data$lat)
data$long <- as.numeric(data$long)

library(letsR)
PAM <- lets.presab.points(data[,1:2], data[,3], xmn=min(data$long), xmx=max(data$long),
                          ymn=min(data$lat), ymx=max(data$lat))
plot(PAM)
plot(PAM, name="Mountain Ringlet")
summary(PAM)
centroids_tp2_rec <- lets.midpoint(PAM, planar=FALSE, method="MCC")
colnames(centroids_tp2_rec) <- c("species","long","lat")
centroids_tp2_rec$cent <- "yes"
data$cent <- "no"

final_centroid <- rbind(data,centroids_tp2_rec)
species <- unique(final_centroid$species)
## plot map for each species with red dot signifying centroid
## colour=cent 
## save graphs to look at
worldmap = map_data('world')

for(i in species) {
  print(i)
  # Printing ggplot within for-loop
  
  temp_plot <- ggplot() + 
    geom_polygon(data = worldmap, 
                 aes(x = long, y = lat, group = group), 
                 fill = 'gray90', color = 'black') + 
    coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
    theme_void() + 
    geom_point(data = final_centroid[final_centroid$species==i,], 
               aes(x = as.numeric(long), 
                   y = as.numeric(lat), colour=cent), size=2) +
    scale_color_manual(values=c("black", "red")) +
    labs(colour = "Centroid value") 
  ggtitle(i) +
    theme(title = element_text(size = 12))
  
  ggsave(temp_plot, file=paste0("../../Maps/Centroids/BNM_centroid_TP2_", i,".png"), width = 20, height = 25, units = "cm")
  Sys.sleep(2)
}

## remove cent and species column
species <- centroids_tp1_rec$species
centroids_tp1_rec <- centroids_tp1_rec[,2:3]
centroids_tp2_rec <- centroids_tp2_rec[,2:3]


## measure distance? and direction of shift
library(geosphere)
x <- bearing(centroids_tp1_rec, centroids_tp2_rec)
## add in time period
centroids_tp1_rec$time_period <- "TP1"
centroids_tp2_rec$time_period <- "TP2"

centroids_rec_hecs <- rbind(centroids_tp1_rec, centroids_tp2_rec)
centroids_rec_hecs$species <- species

direction <- data.frame(bearing=x, species=species)
course <- (x + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
direction$course  <- course
## reorder columns
direction <- direction[,c(2,1,3)]

## plot each species centroid in TP1 and TP2 to check direction
## save graphs to look at
worldmap = map_data('world')

for(i in species) {
  print(i)
  # Printing ggplot within for-loop
  
  temp_plot <- ggplot() + 
    geom_polygon(data = worldmap, 
                 aes(x = long, y = lat, group = group), 
                 fill = 'gray90', color = 'black') + 
    coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
    theme_void() + 
    geom_point(data = centroids_rec_hecs[centroids_rec_hecs$species==i,], 
               aes(x = as.numeric(long), 
                   y = as.numeric(lat), colour=time_period), size=2) +
    geom_line() +
    labs(colour = "Time period") 
  ggtitle(i) +
    theme(title = element_text(size = 12))
  
  ggsave(temp_plot, file=paste0("../../Maps/Centroids/BNM_centroid_shift_", i,".png"), width = 20, height = 25, units = "cm")
  Sys.sleep(2)
}

## get each species classed into 8-point compass directions
## N = 337.5 - 22.5
## NE = 22.5 - 67.5
## E = 67.5 - 112.5
## SE = 112.5 - 157.5
## S = 157.5 - 202.5
## SW = 202.5 - 247.5
## W = 247.5 - 292.5
## NW = 292.5 - 337.5

## long way of doing it but it works
direction$compass_direction <- NULL
NE <- direction %>% filter(between(course,22.5,67.5) )
NE$compass_direction <- "NE"
E <- direction %>% filter(between(course,67.5,112.5) )
E$compass_direction <- "E"
SE <- direction %>% filter(between(course,112.5,157.5) )
SE$compass_direction <- "SE"
S <- direction %>% filter(between(course,157.5,202.5) )
S$compass_direction <- "S"
SW <- direction %>% filter(between(course,202.5,247.5) )
SW$compass_direction <- "SW"
W <- direction %>% filter(between(course,247.5,292.5) )
W$compass_direction <- "W"
NW <- direction %>% filter(between(course,292.5,337.5) )
NW$compass_direction <- "NW"
N <- direction %>% filter(course>337.5 | course<22.5) 
N$compass_direction <- "N"

directions_rec_hecs <- rbind(NE,E,SE,S,SW,W,NW,N)
## nothing has moving in east direction


###############################################################
## well & heavily recorded hectads
rm(list = ls())

bnmdata_well_heavy_hecs <- readRDS("../../Data/BNM/BNMdata_well_heavy_hecs_75_91_12_16.rds")


## Filter 1: northern butterflies only: mountain ringlet, scotch argus, northern brown argus and large heath
northern_butterflies <- c("Large Heath", "Scotch Argus", "Mountain Ringlet", "Northern Brown Argus")
bnmdata_well_heavy_hecs <- bnmdata_well_heavy_hecs[which(bnmdata_well_heavy_hecs$Species %in% northern_butterflies), ]
length(unique(bnmdata_well_heavy_hecs$Species)) ## 4 species

####### Filter 2: species must have enough coverage to calculate range shifts
####### species need to occupy at least 20 hectads in each time period
bnmdata_well_heavy_hecs$Time_period <- ifelse(bnmdata_well_heavy_hecs$Year==2012 | bnmdata_well_heavy_hecs$Year==2013 | bnmdata_well_heavy_hecs$Year==2014
                                       | bnmdata_well_heavy_hecs$Year==2015 | bnmdata_well_heavy_hecs$Year==2016, "TP2", "TP1")
# get rid of the year information, keeping only which hectads each species was recorded in:
bnmdata_well_heavy_hecs2 <- ddply(bnmdata_well_heavy_hecs, .(Hectad,Species,easting,northing,Time_period), summarise,
                                   PRESENCE = 1) 
# and now count the number of hectads per species for each year
tp_species_hecs <- ddply(bnmdata_well_heavy_hecs2, .(Species, Time_period), summarise,
                         HECTADS = sum(PRESENCE))
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_hecs <- tp_species_hecs %>% group_by(Species) %>% filter(all(HECTADS>=20))
species_keep <- unique(tp_species_hecs$Species) ## 3 species (lose mountain ringlet)
## mostly removes rare/restricted species
bnmdata_well_heavy_hecs <-  bnmdata_well_heavy_hecs[which(bnmdata_well_heavy_hecs$Species %in% species_keep),] ## 3 species
## remove duplicated records at 10km resolution - same species/year/hectad (probably have different 1km sites)
bnmdata_well_heavy_hecs <- unique(bnmdata_well_heavy_hecs) ## 1180 rows
length(unique(bnmdata_well_heavy_hecs$Hectad)) ## 288 hectads

## map
worldmap = map_data('world')
well_heavy_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =bnmdata_well_heavy_hecs, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
well_heavy_hecs
## save graph
## scottish coverage very poor - must be lots of hectads in scotland recorded recently (i.e. TP2), but not in TP1
ggsave(well_heavy_hecs, file="../../Maps/Nothern_butterflies_3spp_well_heavy_hecs.png")

### now calculate centroid elevation shift between time periods
## centroid = average elevation of all sampled occurrences during each time period
## and measure the low elevation shift - mean elevation of 10 lowest elevation hectads for each time period

## mean elevation for each species and each time period + SE + confidence intervals
library(Rmisc)
mean_elev_wh_hecs <- summarySE(bnmdata_well_heavy_hecs, measurevar="elevation10x10km", groupvars=c("Species","Time_period"))

## check for normality 
# compute the difference
d <- with(mean_elev_wh_hecs, 
          elevation10x10km[Time_period == "TP1"] - elevation10x10km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.5192 - distribution of the differences are not significantly different from normal distribution 
## i.e. they follow a normal distribution so t-test can be used
qqnorm(d)
qqline(d) ## looks fine

mean_elev_wh <- t.test(elevation10x10km ~ Time_period, data = mean_elev_wh_hecs, paired = TRUE)
mean_elev_wh ## non-significant

## plot anyway
ggpaired(mean_elev_wh_hecs, x = "Time_period", y = "elevation10x10km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Species")+
  xlab("Time period")+
  ylab("Mean elevation")+
  stat_compare_means(method="t.test", paired = TRUE)
## all 3 moved downhill

## create bar plots on elevation change over time 
ggplot(mean_elev_wh_hecs, aes(x=Species, y=elevation10x10km, fill=Time_period)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=elevation10x10km-ci, ymax=elevation10x10km+ci), width=.2, position=position_dodge(0.9))+
  xlab("Common Name") +
  ylab("Mean elevation") +
  theme_classic()

#####################################
## low elevational boundary shift
## this is the 10 lowest OCCUPIED hectads - can't have the same hectad repeated in this 10 just because it was recorded in multiple years

## first remove year info - this means we can only look at occupied hectads within each time period, rather than each year
## e.g. if one hectad is recorded in multiple years within one time period, it would be included multiple times in the 10 lowest hectads
bnmdata_well_heavy_hecs2 <- ddply(bnmdata_well_heavy_hecs, .(Hectad,Species,easting,northing,Time_period,elevation10x10km,lat,lon), summarise,
                           PRESENCE = 1) 
## now order elevation by ascending value
## take the top 10 values 
low_10_elev_wh_hecs <- bnmdata_well_heavy_hecs2 %>% group_by(Species, Time_period) %>% arrange(elevation10x10km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_wh_hecs <- summarySE(low_10_elev_wh_hecs, measurevar="elevation10x10km", groupvars=c("Species","Time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_wh_hecs, 
          elevation10x10km[Time_period == "TP1"] - elevation10x10km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.1457 - distribution is normal
qqnorm(d)
qqline(d) ## looks fine

## paired t-test
low_elev_wh <- t.test(elevation10x10km ~ Time_period, data = low_elev_wh_hecs, paired = TRUE)
low_elev_wh ## significant p=0.03

## plot result
library(ggpubr)
ggpaired(low_elev_wh_hecs, x = "Time_period", y = "elevation10x10km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Species")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="t.test", paired = TRUE)
## all 3 species expand their distribution downhill

## create bar plots on elevation change over time 
ggplot(low_elev_wh_hecs, aes(x=Species, y=elevation10x10km, fill=Time_period)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=elevation10x10km-ci, ymax=elevation10x10km+ci), width=.2, position=position_dodge(0.9))+
  xlab("Common Name") +
  ylab("Low elevation boundary") +
  theme_classic()

########################################
## Latitudinal shifts (trailing edge boundary)
## find 10 most southerly hectads & calculate mean latitude across years

## use same data frame as above that has year info removed
## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_wh_hecs <- bnmdata_well_heavy_hecs2 %>% group_by(Species, Time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_wh_hecs <- summarySE(low_10_lat_wh_hecs, measurevar="lat", groupvars=c("Species","Time_period"))
## higher latitude = more northerly

## check for normality 
# compute the difference
d <- with(low_lat_wh_hecs, 
          lat[Time_period == "TP1"] - lat[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.5932 - normal distribution
qqnorm(d)
qqline(d) ## looks good

## paired t-test
low_lat_wh <- t.test(lat ~ Time_period, data = low_lat_wh_hecs, paired = TRUE)
low_lat_wh ## non-significant

## plot anyway
ggpaired(low_lat_wh_hecs, x = "Time_period", y = "lat",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Species")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="t.test", paired = TRUE)
## not much change

## create bar plots on elevation change over time 
ggplot(low_lat_wh_hecs, aes(x=Species, y=lat, fill=Time_period)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=lat-ci, ymax=lat+ci), width=.2, position=position_dodge(0.9))+
  xlab("Common Name") +
  ylab("Low latititude boundary") +
  theme_classic()

#############################################
## direction of centroid shift 

bnmdata_well_heavy_hecs_tp1 <- bnmdata_well_heavy_hecs[bnmdata_well_heavy_hecs$Time_period=="TP1",]
bnmdata_well_heavy_hecs_tp2 <- bnmdata_well_heavy_hecs[bnmdata_well_heavy_hecs$Time_period=="TP2",]

## first time period
species <- bnmdata_well_heavy_hecs_tp1$Species
long <- bnmdata_well_heavy_hecs_tp1$lon
lat <- bnmdata_well_heavy_hecs_tp1$lat
data <- as.data.frame(cbind(long,lat,species))

data$lat <- as.numeric(data$lat)
data$long <- as.numeric(data$long)

library(letsR)
PAM <- lets.presab.points(data[,1:2], data[,3], xmn=min(data$long), xmx=max(data$long),
                          ymn=min(data$lat), ymx=max(data$lat))
plot(PAM)
plot(PAM, name="Large Heath")
summary(PAM)
centroids_tp1_wh <- lets.midpoint(PAM, planar=FALSE, method="MCC")
colnames(centroids_tp1_wh) <- c("species","long","lat")

## second time period
species <- bnmdata_well_heavy_hecs_tp2$Species
long <- bnmdata_well_heavy_hecs_tp2$lon
lat <- bnmdata_well_heavy_hecs_tp2$lat
data <- as.data.frame(cbind(long,lat,species))

data$lat <- as.numeric(data$lat)
data$long <- as.numeric(data$long)

library(letsR)
PAM <- lets.presab.points(data[,1:2], data[,3], xmn=min(data$long), xmx=max(data$long),
                          ymn=min(data$lat), ymx=max(data$lat))
plot(PAM)
plot(PAM, name="Large Heath")
summary(PAM)
centroids_tp2_wh <- lets.midpoint(PAM, planar=FALSE, method="MCC")
colnames(centroids_tp2_wh) <- c("species","long","lat")

## remove species column (not needed for bearing calculation)
species <- centroids_tp1_wh$species
centroids_tp1_wh <- centroids_tp1_wh[,2:3]
centroids_tp2_wh <- centroids_tp2_wh[,2:3]

## measure direction of shift
library(geosphere)
x <- bearing(centroids_tp1_wh, centroids_tp2_wh)
## create dataframe with bearings
direction <- data.frame(bearing=x, species=species)
course <- (x + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
direction$course  <- course
## reorder columns
direction <- direction[,c(2,1,3)]

## get each species classed into 8-point compass directions
## N = 337.5 - 22.5
## NE = 22.5 - 67.5
## E = 67.5 - 112.5
## SE = 112.5 - 157.5
## S = 157.5 - 202.5
## SW = 202.5 - 247.5
## W = 247.5 - 292.5
## NW = 292.5 - 337.5

## long way of doing it but it works
direction$compass_direction <- NULL
NE <- direction %>% filter(between(course,22.5,67.5) )
NE$compass_direction <- "NE"
E <- direction %>% filter(between(course,67.5,112.5) )
E$compass_direction <- "E"
SE <- direction %>% filter(between(course,112.5,157.5) )
SE$compass_direction <- "SE"
S <- direction %>% filter(between(course,157.5,202.5) )
S$compass_direction <- "S"
SW <- direction %>% filter(between(course,202.5,247.5) )
SW$compass_direction <- "SW"
W <- direction %>% filter(between(course,247.5,292.5) )
W$compass_direction <- "W"
NW <- direction %>% filter(between(course,292.5,337.5) )
NW$compass_direction <- "NW"
N <- direction %>% filter(course>337.5 | course<22.5) 
N$compass_direction <- "N"

directions_well_heavy_hecs <- rbind(NE,E,SE,S,SW,W,NW,N)
## nothing moving in east direction

###############################################################
## heavy hectads
rm(list = ls())

bnmdata_heavy_hecs <- readRDS("../../Data/BNM/BNMdata_heavy_hecs_75_91_12_16.rds")

## Filter 1: northern butterflies only: mountain ringlet, scotch argus, northern brown argus and large heath
northern_butterflies <- c("Large Heath", "Scotch Argus", "Mountain Ringlet", "Northern Brown Argus")
bnmdata_heavy_hecs <- bnmdata_heavy_hecs[which(bnmdata_heavy_hecs$Species %in% northern_butterflies), ]
length(unique(bnmdata_heavy_hecs$Species)) ## 4 species

####### Filter 2: species must have enough coverage to calculate range shifts
####### species need to occupy at least 20 hectads in each time period
bnmdata_heavy_hecs$Time_period <- ifelse(bnmdata_heavy_hecs$Year==2012 | bnmdata_heavy_hecs$Year==2013 | bnmdata_heavy_hecs$Year==2014
                                              | bnmdata_heavy_hecs$Year==2015 | bnmdata_heavy_hecs$Year==2016, "TP2", "TP1")
# get rid of the year information, keeping only which hectads each species was recorded in:
bnmdata_heavy_hecs2 <- ddply(bnmdata_heavy_hecs, .(Hectad,Species,easting,northing,Time_period), summarise,
                                  PRESENCE = 1) 
# and now count the number of hectads per species for each year
tp_species_hecs <- ddply(bnmdata_heavy_hecs2, .(Species, Time_period), summarise,
                         HECTADS = sum(PRESENCE))
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_hecs <- tp_species_hecs %>% group_by(Species) %>% filter(all(HECTADS>=20))
species_keep <- unique(tp_species_hecs$Species) ## 3 species (lose mountain ringlet)
## mostly removes rare/restricted species
bnmdata_heavy_hecs <-  bnmdata_heavy_hecs[which(bnmdata_heavy_hecs$Species %in% species_keep),] ## 3 species
## remove duplicated records at 10km resolution - same species/year/hectad (probably have different 1km sites)
bnmdata_heavy_hecs <- unique(bnmdata_heavy_hecs) ## 734 rows
length(unique(bnmdata_heavy_hecs$Hectad)) ## 134 hectads

## map
worldmap = map_data('world')
heavy_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =bnmdata_heavy_hecs, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
heavy_hecs
## save graph
## scottish coverage very poor - must be lots of hectads in scotland recorded recently (i.e. TP2), but not in TP1
ggsave(heavy_hecs, file="../../Maps/Nothern_butterflies_3spp_heavy_hecs.png")

### now calculate centroid elevation shift between time periods
## centroid = average elevation of all sampled occurrences during each time period
## and measure the low elevation shift - mean elevation of 10 lowest elevation hectads for each time period

## mean elevation for each species and each time period + SE + confidence intervals
library(Rmisc)
mean_elev_heavy_hecs <- summarySE(bnmdata_heavy_hecs, measurevar="elevation10x10km", groupvars=c("Species","Time_period"))

## check for normality 
# compute the difference
d <- with(mean_elev_heavy_hecs, 
          elevation10x10km[Time_period == "TP1"] - elevation10x10km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.7949 - distribution of the differences are not significantly different from normal distribution 
## i.e. they follow a normal distribution so t-test can be used
qqnorm(d)
qqline(d) ## looks fine

mean_elev_heavy <- t.test(elevation10x10km ~ Time_period, data = mean_elev_heavy_hecs, paired = TRUE)
mean_elev_heavy ## non-significant

## plot anyway
ggpaired(mean_elev_heavy_hecs, x = "Time_period", y = "elevation10x10km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Species")+
  xlab("Time period")+
  ylab("Mean elevation")+
  stat_compare_means(method="t.test", paired = TRUE)
## all 3 moved downhill

## create bar plots on elevation change over time 
ggplot(mean_elev_heavy_hecs, aes(x=Species, y=elevation10x10km, fill=Time_period)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=elevation10x10km-ci, ymax=elevation10x10km+ci), width=.2, position=position_dodge(0.9))+
  xlab("Common Name") +
  ylab("Mean elevation") +
  theme_classic()

#####################################
## low elevational boundary shift
## this is the 10 lowest OCCUPIED hectads - can't have the same hectad repeated in this 10 just because it was recorded in multiple years

## first remove year info - this means we can only look at occupied hectads within each time period, rather than each year
## e.g. if one hectad is recorded in multiple years within one time period, it would be included multiple times in the 10 lowest hectads
bnmdata_heavy_hecs2 <- ddply(bnmdata_heavy_hecs, .(Hectad,Species,easting,northing,Time_period,elevation10x10km,lat,lon), summarise,
                                  PRESENCE = 1) 
## now order elevation by ascending value
## take the top 10 values 
low_10_elev_heavy_hecs <- bnmdata_heavy_hecs2 %>% group_by(Species, Time_period) %>% arrange(elevation10x10km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_heavy_hecs <- summarySE(low_10_elev_heavy_hecs, measurevar="elevation10x10km", groupvars=c("Species","Time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_heavy_hecs, 
          elevation10x10km[Time_period == "TP1"] - elevation10x10km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.83 - distribution is normal
qqnorm(d)
qqline(d) ## looks fine

## paired t-test
low_elev_heavy <- t.test(elevation10x10km ~ Time_period, data = low_elev_heavy_hecs, paired = TRUE)
low_elev_heavy ## significant p=0.019

## plot result
ggpaired(low_elev_heavy_hecs, x = "Time_period", y = "elevation10x10km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Species")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="t.test", paired = TRUE)
## all 3 species expand their distribution downhill

## create bar plots on elevation change over time 
ggplot(low_elev_heavy_hecs, aes(x=Species, y=elevation10x10km, fill=Time_period)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=elevation10x10km-ci, ymax=elevation10x10km+ci), width=.2, position=position_dodge(0.9))+
  xlab("Common Name") +
  ylab("Low elevation boundary") +
  theme_classic()

########################################
## Latitudinal shifts (trailing edge boundary)
## find 10 most southerly hectads & calculate mean latitude across years

## use same data frame as above that has year info removed
## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_heavy_hecs <- bnmdata_heavy_hecs2 %>% group_by(Species, Time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_heavy_hecs <- summarySE(low_10_lat_heavy_hecs, measurevar="lat", groupvars=c("Species","Time_period"))
## higher latitude = more northerly

## check for normality 
# compute the difference
d <- with(low_lat_heavy_hecs, 
          lat[Time_period == "TP1"] - lat[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.3129 - normal distribution
qqnorm(d)
qqline(d) ## looks good

## paired t-test
low_lat_heavy <- t.test(lat ~ Time_period, data = low_lat_heavy_hecs, paired = TRUE)
low_lat_heavy ## non-significant

## plot anyway
ggpaired(low_lat_heavy_hecs, x = "Time_period", y = "lat",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Species")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="t.test", paired = TRUE)
## not much change

## create bar plots on elevation change over time 
ggplot(low_lat_heavy_hecs, aes(x=Species, y=lat, fill=Time_period)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=lat-ci, ymax=lat+ci), width=.2, position=position_dodge(0.9))+
  xlab("Common Name") +
  ylab("Low latititude boundary") +
  theme_classic()


#############################################
## direction of centroid shift 

bnmdata_heavy_hecs_tp1 <- bnmdata_heavy_hecs[bnmdata_heavy_hecs$Time_period=="TP1",]
bnmdata_heavy_hecs_tp2 <- bnmdata_heavy_hecs[bnmdata_heavy_hecs$Time_period=="TP2",]

## first time period
species <- bnmdata_heavy_hecs_tp1$Species
long <- bnmdata_heavy_hecs_tp1$lon
lat <- bnmdata_heavy_hecs_tp1$lat
data <- as.data.frame(cbind(long,lat,species))

data$lat <- as.numeric(data$lat)
data$long <- as.numeric(data$long)

library(letsR)
PAM <- lets.presab.points(data[,1:2], data[,3], xmn=min(data$long), xmx=max(data$long),
                          ymn=min(data$lat), ymx=max(data$lat))
plot(PAM)
plot(PAM, name="Large Heath")
summary(PAM)
centroids_tp1_heavy <- lets.midpoint(PAM, planar=FALSE, method="MCC")
colnames(centroids_tp1_heavy) <- c("species","long","lat")

## second time period
species <- bnmdata_heavy_hecs_tp2$Species
long <- bnmdata_heavy_hecs_tp2$lon
lat <- bnmdata_heavy_hecs_tp2$lat
data <- as.data.frame(cbind(long,lat,species))

data$lat <- as.numeric(data$lat)
data$long <- as.numeric(data$long)

library(letsR)
PAM <- lets.presab.points(data[,1:2], data[,3], xmn=min(data$long), xmx=max(data$long),
                          ymn=min(data$lat), ymx=max(data$lat))
plot(PAM)
plot(PAM, name="Large Heath")
summary(PAM)
centroids_tp2_heavy <- lets.midpoint(PAM, planar=FALSE, method="MCC")
colnames(centroids_tp2_heavy) <- c("species","long","lat")

## remove species column (not needed for bearing calculation)
species <- centroids_tp1_heavy$species
centroids_tp1_heavy <- centroids_tp1_heavy[,2:3]
centroids_tp2_heavy <- centroids_tp2_heavy[,2:3]

## measure direction of shift
library(geosphere)
x <- bearing(centroids_tp1_heavy, centroids_tp2_heavy)
## create dataframe with bearings
direction <- data.frame(bearing=x, species=species)
course <- (x + 360) %% 360 # add full circle, i.e. +360, and determine modulo for 360
direction$course  <- course
## reorder columns
direction <- direction[,c(2,1,3)]

## get each species classed into 8-point compass directions
## N = 337.5 - 22.5
## NE = 22.5 - 67.5
## E = 67.5 - 112.5
## SE = 112.5 - 157.5
## S = 157.5 - 202.5
## SW = 202.5 - 247.5
## W = 247.5 - 292.5
## NW = 292.5 - 337.5

## long way of doing it but it works
direction$compass_direction <- NULL
NE <- direction %>% filter(between(course,22.5,67.5) )
NE$compass_direction <- "NE"
E <- direction %>% filter(between(course,67.5,112.5) )
E$compass_direction <- "E"
SE <- direction %>% filter(between(course,112.5,157.5) )
SE$compass_direction <- "SE"
S <- direction %>% filter(between(course,157.5,202.5) )
S$compass_direction <- "S"
SW <- direction %>% filter(between(course,202.5,247.5) )
SW$compass_direction <- "SW"
W <- direction %>% filter(between(course,247.5,292.5) )
W$compass_direction <- "W"
NW <- direction %>% filter(between(course,292.5,337.5) )
NW$compass_direction <- "NW"
N <- direction %>% filter(course>337.5 | course<22.5) 
N$compass_direction <- "N"

directionsheavy_hecs <- rbind(NE,E,SE,S,SW,W,NW,N)
## nothing moving in east direction






# 
# 
# 
# ## put results in one dataframe 
# recorded_results <- merge(mean_elev_rec, low_elev_rec, by="Species", all=T)
# recorded_results <- merge(recorded_results, low_lat_rec, by="Species", all=T)
# recorded_results$recording_level <- "Recorded"
# 
# well_rec_results <- merge(mean_elev_well, low_elev_well, by="Species", all=T)
# well_rec_results <- merge(well_rec_results, low_lat_well, by="Species", all=T)
# well_rec_results$recording_level <- "Well Recorded"
# 
# heavy_rec_results <- merge(mean_elev_heavy, low_elev_heavy, by="Species", all=T)
# heavy_rec_results <- merge(heavy_rec_results, low_lat_heavy, by="Species", all=T)
# heavy_rec_results$recording_level <- "Heavily Recorded"
# 
# bnm_range_shift_results <- rbind(recorded_results, well_rec_results, heavy_rec_results)
# ## change to long format
# bnm_range_shift_results2 <- reshape(bnm_range_shift_results, direction="long", 
#         varying=list(c("mean_elev_tp1","low_elev_tp1", "low_lat_tp1"), c("mean_elev_tp2","low_elev_tp2", "low_lat_tp2"),
#                      c("mean_elev_diff","low_elev_diff", "low_lat_diff")), 
#         v.names=c("TP1","TP2","Difference"))
# bnm_range_shift_results2$time[bnm_range_shift_results2$time==1] <- "Mean elevation"
# bnm_range_shift_results2$time[bnm_range_shift_results2$time==2] <- "Low elevation"
# bnm_range_shift_results2$time[bnm_range_shift_results2$time==3] <- "Low latitude"
# bnm_range_shift_results2$id <- NULL
# write.csv(bnm_range_shift_results2, file="../../Results/bnm_range_shift_results.csv", row.names=FALSE)
# 
