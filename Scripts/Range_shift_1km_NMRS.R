##########################
#### user: Lisbeth Hordley
#### date: September 2021
#### info: Range shift analysis at 1km scale using NMRS data

rm(list = ls())

## load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(sf)
library(ggplot2)
library(plyr)

## read in nmrs data
nmrsdata <- fread("../../Range shift/Data/NMRS/NMRS_Atlas_Trend_Data_20190311_final.csv", stringsAsFactors = FALSE)

## tidy data
# Obtain year from date
nmrsdata$Date <- dmy(nmrsdata$Date)
nmrsdata$Year <- year(nmrsdata$Date)
# Remove unecessary columns
nmrsdata[,c("Key","Code","Abundance","Life Stage","Date","Vice-county","Easting","Northing","_10k"):=NULL] 
# remove easting and northing because they are re-calculated below based on the 1km gridref
## this is just to make sure they are all at the same scale (i.e. taken from the bottom left hand corner of square)
head(nmrsdata)
colnames(nmrsdata)[8] <- "Hectad"
nmrsdata <- nmrsdata %>% mutate_all(na_if,"") ## change blanks to NAs
colSums(is.na(nmrsdata))
x <- nmrsdata[is.na(nmrsdata$Grid_square),]
x_tp1 <- x[x$Year>=1975&x$Year<=1991,]
x_tp2 <- x[x$Year>=2012&x$Year<=2016,]

## 45079 NAs for gridsquares - these records do not have 1km resolution 
## remove these records
nmrsdata <- na.omit(nmrsdata) ## 20533692 rows

## convert 1km gridrefs to lat and lon and Easting and Northing
library(rnrfa)
gridrefs <- as.vector(unique(nmrsdata[,3])) ## 62539 1km gridrefs
lon_lat <- as.data.frame(osg_parse(grid_refs=gridrefs$Grid_square, coord_system = "WGS84"))
east_north <- as.data.frame(osg_parse(grid_refs=gridrefs$Grid_square))
gridref_latlon_eastnorth <- cbind(gridrefs, lon_lat, east_north)
## convert gridref to east/north - gives easting and northing at 10km scale - used to measure distances between sites later
## merge back into nmrs data
nmrsdata <- merge(nmrsdata, gridref_latlon_eastnorth, by="Grid_square")
head(nmrsdata)
## save NMRS cleaned data
saveRDS(nmrsdata, "../../Range shift/Data/NMRS/NMRS_grid_square_cleaned.rds")

## elevation data
library(raster)
tmp = raster("../../Data/elevation1x1_new.tif")
crs(tmp) ## projection
res(tmp) ## resolution (1000 x 1000m, 1x1km)
plot(tmp)
hist(values(tmp))

## re-project
tmp_newproj <- projectRaster(tmp,crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(tmp_newproj, axes = T)
tmp_elev <- raster::extract(tmp_newproj, gridref_latlon_eastnorth[,2:3], df=T)
tmp_elev <- cbind(tmp_elev, gridref_latlon_eastnorth)
## multiply elevation by 10 - see metadata https://www.eea.europa.eu/data-and-maps/data/digital-elevation-model-of-europe?fbclid=IwAR1yqtoCIjiV1QIe30h8DCaYhWH-7GO0X68CKLAkD5FY9e5WMs6OPU9YSas 
tmp_elev$elevation1x1_new <- tmp_elev$elevation1x1_new*10 ## elevation in meters 
tmp_elev <- tmp_elev[,-1]
colnames(tmp_elev)[1] <- "elevation1x1km"

## map of elevation at NMRS 1km sqaures/points
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = tmp_elev, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour = elevation1x1km), size=1) + 
  scale_color_viridis_c(name="Elevation") + 
  theme_void() +
  theme(title = element_text(size = 12))

nmrsdata <- readRDS("../../Range shift/Data/NMRS/NMRS_grid_square_cleaned.rds") ## 20,533,692 rows
## remove lat/lon and easting/northing - not needed to merge in, just use gridref
tmp_elev <- tmp_elev[,-c(3:6)]
## if throwing memory issues: memory.limit(size = 15000)
nmrsdata <- merge(nmrsdata, tmp_elev, by="Grid_square", all=T)
## 20,533,692 rows in NMRS data (same as when first read in)
head(nmrsdata)
## save data
saveRDS(nmrsdata, file="../../Range shift/Data/NMRS/NMRS_grid_square_elevation.rds")


#############################################################################################################
############## look at data in more detail
rm(list = ls())
options(scipen=999)

nmrsdata <- readRDS("../../Range shift/Data/NMRS/NMRS_grid_square_elevation.rds")
## easting/northing, lat/lon and elevation are at 1km resolution

## TP2 will be most recent 5 year period (2012-2016)
## TP1 needs to be at least 10 or 15? years before
## and be a 5 year period with the most grid squares in total

## try a 10 year gap first (1998-2002 and earlier)
## do a 5 year rolling window from 1990 to 1998 and see which has highest number of grid squares 

year.list <- 1990:1998
final <- NULL
for (i in year.list){ # loop through years
  start.year<-i
  end.year<-i+4
  nmrs_5_year_data<-nmrsdata[nmrsdata$Year>=start.year&nmrsdata$Year<=end.year,]
  grid_squares <- unique(nmrs_5_year_data$Grid_square)
  tmp <- data.frame(no_grid_squares=length(grid_squares), start_year=i, end_year=end.year)
  final <- rbind(tmp, final)
  }
## closest (1998-2002) 5 year period has the highest number of grid squares (12777)
## 9567 grid squares if we left a 15 year gap
nmrs_recent<-nmrsdata[nmrsdata$Year>=2012&nmrsdata$Year<=2016,]
length(unique(nmrs_recent$Grid_square)) ## 31254 grid squares in most recent 5 year period

# count the number of recorded hectads and the number of species*hectad records in each year 
library(plyr)
nmrsdata$PRESENCE <- 1
nmrs_rec_gs_years <- ddply(nmrsdata, .(Grid_square, Year, easting, northing), numcolwise(mean)) ## removes species data, just grid squares recorded for each year

### time periods:
tp1 <- 1998:2002
tp2 <- 2012:2016

# first follow Hickling et al directly by calculating which grid squares were recorded at least once in each interval
# first select out the intervals
nmrs_rec_gs_years_tp1 <- nmrs_rec_gs_years[which(nmrs_rec_gs_years$Year %in% tp1), ]
nmrs_rec_gs_years_tp2 <- nmrs_rec_gs_years[which(nmrs_rec_gs_years$Year %in% tp2), ]

# then get rid of the year information to leave a non-redundant list of grid squares recorded in each interval
nmrs_rec_gs_tp1 <- ddply(nmrs_rec_gs_years_tp1, .(Grid_square), numcolwise(mean))
nmrs_rec_gs_tp2 <- ddply(nmrs_rec_gs_years_tp2, .(Grid_square), numcolwise(mean))
## this just takes the mean year for each grid square - just gives a unique list of grid squares recorded in each time period
## 21,270 in tp1
## 55,232 in tp2

# plot these on a map
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrs_rec_gs_tp2,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat)), size=1, colour="blue") +
  geom_point(data =nmrs_rec_gs_tp1,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat)), size=1, colour="red") +
  theme(title = element_text(size = 12))
## much better coverage of Scotland in TP2

## now we want to recombine these datasets in a way that ditches grid squares not present in both TPs
nmrs_rec_gs_tp1_vec <- as.character(nmrs_rec_gs_tp1$Grid_square)
nmrs_rec_gs_tp2_vec <- as.character(nmrs_rec_gs_tp2$Grid_square)

nmrs_rec_gs_full <- append(nmrs_rec_gs_tp1_vec,nmrs_rec_gs_tp2_vec)

nmrs_rec_gs_full <- data.frame(cbind(nmrs_rec_gs_full,1)) ## list of all grid squares in both TPs
colnames(nmrs_rec_gs_full) <- c("Grid_square","Intervals")
nmrs_rec_gs_full$Grid_square <- as.factor(as.character(nmrs_rec_gs_full$Grid_square)) ## 44031
nmrs_rec_gs_full$Intervals <- as.numeric(as.character(nmrs_rec_gs_full$Intervals))


nmrs_rec_gs_good <- ddply(nmrs_rec_gs_full, .(Grid_square), numcolwise(sum)) ## grid squares either have 1 (only recorded in one TP), or 2 (recorded in both TPs)
## 38048 unique hectads
nmrs_rec_gs_good <- nmrs_rec_gs_good[which(nmrs_rec_gs_good$Intervals > 1), ] ## only keep grid squares recorded in both (i.e. intervals > 1)
## 5983 good grid squares

# drop unused levels
nmrs_rec_gs_good <- droplevels(nmrs_rec_gs_good)

# generate the vector
nmrs_good_gs <- levels(nmrs_rec_gs_good$Grid_square)

# pick out the good records
length(unique(nmrsdata$Grid_square)) ## 62539 total grid squares
nmrsdata_good_gs <- nmrsdata[which(nmrsdata$Grid_square %in% nmrs_good_gs), ]
length(unique(nmrsdata_good_gs$Grid_square)) ## 5983
# we lose over 90% of grid squares...
# or 85% of unique grid squares recorded in years of interest 

nmrs_good_gs2 <- ddply(nmrsdata_good_gs, .(Grid_square, Year), numcolwise(mean))

plot <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrs_good_gs2, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), size=1, colour="purple") +
  theme(title = element_text(size = 12))
plot

## now further restrict based on recording effort
# start with the list of recorded grid squares, re-attach the eastings and northings
nmrs_good <- ddply(nmrs_rec_gs_years, .(Grid_square,easting,northing), summarise,
                   COUNT = sum(PRESENCE))

# select out only the good hectads
nmrs_good <- nmrs_good[which(nmrs_good$Grid_square %in% nmrs_good_gs), ]
time_period1 <- as.data.frame(tp1)
colnames(time_period1)[1] <- "year"
time_period1$interval <- "tp1"
time_period2 <- as.data.frame(tp2)
colnames(time_period2)[1] <- "year"
time_period2$interval <- "tp2"
time_period <- rbind(time_period1, time_period2)

nmrs_good$easting <- as.numeric(nmrs_good$easting)
nmrs_good$northing <- as.numeric(nmrs_good$northing)

# now start to construct the loop:
gs_records <- data.frame(Year = numeric(),
                          Grid_square = factor(),
                          Easting = numeric(),
                          Northing = numeric(),
                          gs_SR = numeric(),
                          reg_SR = numeric(),
                          gs_perc = numeric())
i <- 1

for (x in nmrs_good_gs){
  if (round(i, -1) == i){
    print(i)
  }
  gs <- nmrs_good[which(nmrs_good$Grid_square == x), ]    # pick out details of focal hectad
  candidates <- nmrs_good[which(nmrs_good$Grid_square != x), ]   # pick out details of all others
  
  gs_east <- gs[[1,2]]     # easting of target
  gs_north <- gs[[1,3]]    # northing of target
  
  candidates$east_diff <- candidates$easting - gs_east             # longitudinal difference
  candidates$north_diff <- candidates$northing - gs_north          # latitudinal difference
  
  candidates$distance <- sqrt((candidates$east_diff^2) + (candidates$north_diff^2))    # absolute difference
  
  candidates <- candidates[order(candidates$distance),] # sort by distance ascending
  
  closest <- candidates[1:100,] # select out closest 100
  
  # now we want to calculate species richness from the hectad and from the region in both time periods
  for (n in unique(time_period$interval)){
    y <- time_period$year[time_period$interval==n] ## tp years of interest
    year_recs <- nmrsdata_good_gs[nmrsdata_good_gs$Year %in% y, ] ## nmrs good hectad data for tp of interest
    
    gs_recs <- year_recs[which(year_recs$Grid_square == x), ] # pull out hectad records
    gs_SR <- length(unique(gs_recs$Common_Name))       # calculate species richness
    
    reg_recs <- year_recs[which(year_recs$Grid_square %in% closest$Grid_square), ]  # pull out region records
    reg_recs <- rbind(reg_recs,gs_recs)    # add in hectad records (they're part of the regional richness too!)
    reg_SR <- length(unique(reg_recs$Common_Name))
    
    gs_perc <- gs_SR*100/reg_SR
    
    out <- cbind(x,n,gs_east,gs_north,gs_SR,reg_SR,gs_perc)
    gs_records <- rbind(gs_records,out)
  }    
  i <- i+1
}

colnames(gs_records) <- c("GRID_SQUARE","TIME_PERIOD","EASTING","NORTHING","GRID_SQUARE_SR","REGION_SR","PERCENT_RECORDED")
summary(gs_records)

# make things that should be numeric, numeric
gs_records$EASTING <- as.numeric(as.character(gs_records$EASTING))
gs_records$NORTHING <- as.numeric(as.character(gs_records$NORTHING))
gs_records$GRID_SQUARE_SR <- as.numeric(as.character(gs_records$GRID_SQUARE_SR))
gs_records$REGION_SR <- as.numeric(as.character(gs_records$REGION_SR))
gs_records$PERCENT_RECORDED <- as.numeric(as.character(gs_records$PERCENT_RECORDED))

summary(gs_records)

# most of the hectad*year combinations are "well recorded" but only slightly under a quarter are "heavily recorded"

# label each hectad according to its recording level in each time period
gs_records$RECORDING.LEVEL <- as.factor(ifelse(gs_records$PERCENT_RECORDED >= 25, "Heavily recorded",
                                                ifelse(gs_records$PERCENT_RECORDED >= 10, "Well recorded",
                                                       ifelse(gs_records$PERCENT_RECORDED >0, "Recorded",
                                                              "Not recorded"))))
## get lat/lon values again to plot with ggplot
lat_lon <- nmrsdata_good_gs %>% distinct(Grid_square, lat, lon, .keep_all = FALSE)

gs_records <- merge(gs_records, lat_lon, by.x="GRID_SQUARE", by.y="Grid_square")

# this data takes a while to generate so let's back it up
write.csv(gs_records, "../../Data/NMRS/Grid_square_recording_levels_1975_1991_2012_2016.csv", row.names = F)
gs_records <- read.csv("../../Data/NMRS/Grid_square_recording_levels_1975_1991_2012_2016.csv", header=TRUE)

## Analysis 1: hectads with at least one species recorded in BOTH time periods: 5983 grid squares
## file = nmrsdata_good_gs (needs subsetting by years)
years <- c(1998:2002, 2012:2016)
nmrsdata_rec_gs <- nmrsdata_good_gs[which(nmrsdata_good_gs$Year %in% years), ]
length(unique(nmrsdata_rec_gs$Grid_square)) ## 5983 hectads
length(unique(nmrsdata_rec_gs$Year)) ## 10 years
## save file
saveRDS(nmrsdata_rec_gs, "../../Data/NMRS/NMRSdata_rec_gs_98_02_12_16.rds")

## Analysis 2a: hectads classed as well-recorded in BOTH time periods: 143 GRID SQUARES (not doable)
well_gs <- gs_records %>%
  filter(RECORDING.LEVEL == "Well recorded") %>%
  group_by(GRID_SQUARE) %>%
  dplyr::summarise(n_row=n())
well_gs <- well_gs[which(well_gs$n_row > 1), ] 
well_gs$n_row <- NULL ## 143

## Analysis 2b: hectads classed as heavily recorded in BOTH time periods: 303 grid squares
heavy_gs <- gs_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded") %>%
  group_by(GRID_SQUARE) %>%
  dplyr::summarise(n_row=n())
heavy_gs <- heavy_gs[which(heavy_gs$n_row > 1), ] 
heavy_gs$n_row <- NULL ## 303 grid squares
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_heavy_gs <- nmrsdata[which(nmrsdata$Grid_square %in% heavy_gs$GRID_SQUARE)]
length(unique(nmrsdata_heavy_gs$Grid_square)) ## 303 grid squares
nmrsdata_heavy_gs <- nmrsdata_heavy_gs[which(nmrsdata_heavy_gs$Year %in% years)]
length(unique(nmrsdata_heavy_gs$Year)) ## 10 years
## save file
saveRDS(nmrsdata_heavy_gs, "../../Data/NMRS/NMRSdata_heavy_gs_98_02_12_16.rds")

## Analysis 2c: hectads classed as well OR heavily recorded in BOTH time periods: 732 grid squares
well_heavy_gs <- gs_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded" | RECORDING.LEVEL=="Well recorded") %>%
  group_by(GRID_SQUARE) %>%
  dplyr::summarise(n_row=n())
well_heavy_gs <- well_heavy_gs[which(well_heavy_gs$n_row > 1), ] 
well_heavy_gs$n_row <- NULL
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_wh_gs <- nmrsdata[which(nmrsdata$Grid_square %in% well_heavy_gs$GRID_SQUARE)]
length(unique(nmrsdata_wh_gs$Grid_square)) ## 732 grid squares
nmrsdata_wh_gs <- nmrsdata_wh_gs[which(nmrsdata_wh_gs$Year %in% years)]
length(unique(nmrsdata_wh_gs$Year)) ## 10 years
## save file
saveRDS(nmrsdata_wh_gs, "../../Data/NMRS/NMRSdata_well_heavy_gs_98_02_12_16.rds")


########################################################################################
## for now do analysis on these 5983 grid squares with no recording effort accounted for

## recorded hectads
rm(list = ls())

nmrsdata_rec_gs <- readRDS("../../Data/NMRS/NMRSdata_rec_gs_98_02_12_16.rds")

## Filter 1: upland moth species only (either obligate or predominate use of uplands)
trait_data <- read.csv("../../Data/ecological_traits.csv", header=TRUE)

length(unique(nmrsdata_rec_gs$Common_Name)) ## 823 species 

## subset species by upland levels 1 and 2 
north_species <- subset(trait_data, X7_montane_upland == 1 | X7_montane_upland == 2)
## remove nymphalids (butterflies)
north_species <- north_species[!(north_species$family=="Nymphalidae"),] ## 36 species
north_species <- north_species[,c(3,4,179)]
#write.csv(north_species, "../../Data/north_species2.csv", row.names=F)

## subset species in nmrs data
nmrsdata_rec_gs <- merge(nmrsdata_rec_gs, north_species, by.x="Preferred_Name", by.y="scientific_name")
length(unique(nmrsdata_rec_gs$Common_Name)) ## 34 species

####### Filter 2: species must have enough coverage to calculate range shifts
####### Option 2: species need to occupy at least 20 hectads in each time period
nmrsdata_rec_gs$Time_period <- ifelse(nmrsdata_rec_gs$Year==2012 | nmrsdata_rec_gs$Year==2013 | nmrsdata_rec_gs$Year==2014
                                        | nmrsdata_rec_gs$Year==2015 | nmrsdata_rec_gs$Year==2016, "TP2", "TP1")
# get rid of the year information, keeping only which hectads each species was recorded in:
nmrsdata_rec_gs2 <- ddply(nmrsdata_rec_gs, .(Grid_square,Preferred_Name,Common_Name,easting,northing,Time_period), summarise,
                            PRESENCE = 1) 
# and now count the number of hectads per species for each year
tp_species_gs <- ddply(nmrsdata_rec_gs2, .(Preferred_Name,Common_Name, Time_period), summarise,
                         Grid_squares = sum(PRESENCE))
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_gs <- tp_species_gs %>% group_by(Common_Name) %>% filter(all(Grid_squares>=20))
species_keep <- unique(tp_species_gs$Common_Name) ## 14 species 
## mostly removes rare/restricted species
nmrsdata_rec_gs <-  nmrsdata_rec_gs[which(nmrsdata_rec_gs$Common_Name %in% species_keep),] ## 14 species
## remove duplicated records at 10km resolution - same species/year/hectad (probably have different 1km sites)
nmrsdata_rec_gs <- unique(nmrsdata_rec_gs) ## 2893
length(unique(nmrsdata_rec_gs$Grid_square)) ## 793 grid squares

## map
worldmap = map_data('world')
rec_gs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_rec_gs, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
rec_gs
ggsave(rec_gs, file="../../Maps/Upland_moth_14spp_rec_gs.png")

### now calculate centroid elevation shift between time periods
## centroid = average elevation of all sampled occurrences during each time period
## and measure the low elevation shift - mean elevation of 10 lowest elevation hectads for each time period

## mean elevation for each species and each time period + SE + confidence intervals
library(Rmisc)
mean_elev_rec_gs <- summarySE(nmrsdata_rec_gs, measurevar="elevation1x1km", groupvars=c("Common_Name","Time_period"))

## check for normality 
# compute the difference
d <- with(mean_elev_rec_gs, 
          elevation1x1km[Time_period == "TP1"] - elevation1x1km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.3069 - distribution of the differences are not significantly different from normal distribution 
## i.e. they follow a normal distribution so t-test can be used
qqnorm(d)
qqline(d) ## looks fine
## northern dart, Rannoch Looper, Weaver's Wave
# mean_elev_rec_gs <- mean_elev_rec_gs[!(mean_elev_rec_gs$Common_Name=="Northern Dart"),]
# mean_elev_rec_gs <- mean_elev_rec_gs[!(mean_elev_rec_gs$Common_Name=="Rannoch Looper"),]
# mean_elev_rec_gs <- mean_elev_rec_gs[!(mean_elev_rec_gs$Common_Name=="Weaver's Wave"),]

mean_elev_rec <- t.test(elevation1x1km ~ Time_period, data = mean_elev_rec_gs, paired = TRUE)
mean_elev_rec ## non-significant

## plot anyway
library(ggpubr)
ggpaired(mean_elev_rec_gs, x = "Time_period", y = "elevation1x1km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Mean elevation")+
  stat_compare_means(method="t.test", paired = TRUE)
## non-significant downwards shift..

#####################################
## low elevational boundary shift
## this is the 10 lowest OCCUPIED hectads - can't have the same hectad repeated in this 10 just because it was recorded in multiple years

## first remove year info - this means we can only look at occupied hectads within each time period, rather than each year
## e.g. if one hectad is recorded in multiple years within one time period, it would be included multiple times in the 10 lowest hectads
nmrsdata_rec_gs2 <- ddply(nmrsdata_rec_gs, .(Grid_square,Common_Name,easting,northing,Time_period,elevation1x1km,lat,lon), summarise,
                            PRESENCE = 1) 
## now order elevation by ascending value
## and take the top 10 values 
low_10_elev_rec_gs <- nmrsdata_rec_gs2 %>% group_by(Common_Name, Time_period) %>% arrange(elevation1x1km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_rec_gs <- summarySE(low_10_elev_rec_gs, measurevar="elevation1x1km", groupvars=c("Common_Name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_rec_gs, 
          elevation1x1km[Time_period == "TP1"] - elevation1x1km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.8782 - distribution is normal
qqnorm(d)
qqline(d) ## good

## wilcoxon signed-rank test instead
low_elev_rec <- t.test(elevation1x1km ~ Time_period, data = low_elev_rec_gs, paired = TRUE)
low_elev_rec ## significant

## plot result
ggpaired(low_elev_rec_gs, x = "Time_period", y = "elevation1x1km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="t.test", paired = TRUE)
## low elevation boundary moving DOWN
## golden rod brindle isn't in this analysis

########################################
## Latitudinal shifts (trailing edge boundary)

## find 10 most southerly hectads & calculate mean latitude across years

## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_rec_gs <- nmrsdata_rec_gs2 %>% group_by(Common_Name, Time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_rec_gs <- summarySE(low_10_lat_rec_gs, measurevar="lat", groupvars=c("Common_Name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_lat_rec_gs, 
          lat[Time_period == "TP1"] - lat[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.07811 - normal distribution
qqnorm(d)
qqline(d) ## looks good

## wilcoxon signed-rank test instead
low_lat_rec <- t.test(lat ~ Time_period, data = low_lat_rec_gs, paired = TRUE)
low_lat_rec ## non-significant

## plot anyway
ggpaired(low_lat_rec_gs, x = "Time_period", y = "lat",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="t.test", paired = TRUE)
## very little change over time



## well & heavily recorded grid squares
rm(list = ls())

nmrsdata_wh_gs <- readRDS("../../Data/NMRS/NMRSdata_well_heavy_gs_98_02_12_16.rds")

## Filter 1: upland moth species only (either obligate or predominate use of uplands)
trait_data <- read.csv("../../Data/ecological_traits.csv", header=TRUE)

length(unique(nmrsdata_rec_gs$Common_Name)) ## 823 species 

## subset species by upland levels 1 and 2 
north_species <- subset(trait_data, X7_montane_upland == 1 | X7_montane_upland == 2)
## remove nymphalids (butterflies)
north_species <- north_species[!(north_species$family=="Nymphalidae"),] ## 36 species
north_species <- north_species[,c(3,4,179)]
#write.csv(north_species, "../../Data/north_species2.csv", row.names=F)

## subset species in nmrs data
nmrsdata_wh_gs <- merge(nmrsdata_wh_gs, north_species, by.x="Preferred_Name", by.y="scientific_name")
length(unique(nmrsdata_wh_gs$Common_Name)) ## 26 species

####### Filter 2: species must have enough coverage to calculate range shifts
####### Option 2: species need to occupy at least 20 hectads in each time period
nmrsdata_wh_gs$Time_period <- ifelse(nmrsdata_wh_gs$Year==2012 | nmrsdata_wh_gs$Year==2013 | nmrsdata_wh_gs$Year==2014
                                      | nmrsdata_wh_gs$Year==2015 | nmrsdata_wh_gs$Year==2016, "TP2", "TP1")
# get rid of the year information, keeping only which hectads each species was recorded in:
nmrsdata_wh_gs2 <- ddply(nmrsdata_wh_gs, .(Grid_square,Preferred_Name,Common_Name,easting,northing,Time_period), summarise,
                          PRESENCE = 1) 
# and now count the number of hectads per species for each year
tp_species_gs <- ddply(nmrsdata_wh_gs2, .(Preferred_Name,Common_Name, Time_period), summarise,
                       Grid_squares = sum(PRESENCE))
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_gs <- tp_species_gs %>% group_by(Common_Name) %>% filter(all(Grid_squares>=20))
species_keep <- unique(tp_species_gs$Common_Name) ## 5 species 
## mostly removes rare/restricted species
nmrsdata_wh_gs <-  nmrsdata_wh_gs[which(nmrsdata_wh_gs$Common_Name %in% species_keep),] ## 5 species
nmrsdata_wh_gs <- unique(nmrsdata_wh_gs) ## 745
length(unique(nmrsdata_wh_gs$Grid_square)) ## 193 grid squares

## map
worldmap = map_data('world')
wh_gs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_wh_gs, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
wh_gs
ggsave(wh_gs, file="../../Maps/Upland_moth_5spp_wh_gs_1km.png")

### now calculate centroid elevation shift between time periods
## centroid = average elevation of all sampled occurrences during each time period
## and measure the low elevation shift - mean elevation of 10 lowest elevation hectads for each time period

## mean elevation for each species and each time period + SE + confidence intervals
library(Rmisc)
mean_elev_wh_gs <- summarySE(nmrsdata_wh_gs, measurevar="elevation1x1km", groupvars=c("Common_Name","Time_period"))

## check for normality 
# compute the difference
d <- with(mean_elev_wh_gs, 
          elevation1x1km[Time_period == "TP1"] - elevation1x1km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.3069 - distribution of the differences are not significantly different from normal distribution 
## i.e. they follow a normal distribution so t-test can be used
qqnorm(d)
qqline(d) ## looks fine
## northern dart, Rannoch Looper, Weaver's Wave
# mean_elev_rec_gs <- mean_elev_rec_gs[!(mean_elev_rec_gs$Common_Name=="Northern Dart"),]
# mean_elev_rec_gs <- mean_elev_rec_gs[!(mean_elev_rec_gs$Common_Name=="Rannoch Looper"),]
# mean_elev_rec_gs <- mean_elev_rec_gs[!(mean_elev_rec_gs$Common_Name=="Weaver's Wave"),]

mean_elev_wh <- t.test(elevation1x1km ~ Time_period, data = mean_elev_wh_gs, paired = TRUE)
mean_elev_wh ## non-significant

## plot anyway
library(ggpubr)
ggpaired(mean_elev_wh_gs, x = "Time_period", y = "elevation1x1km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Mean elevation")+
  stat_compare_means(method="t.test", paired = TRUE)
## non-significant downwards shift..

#####################################
## low elevational boundary shift
## this is the 10 lowest OCCUPIED hectads - can't have the same hectad repeated in this 10 just because it was recorded in multiple years

## first remove year info - this means we can only look at occupied hectads within each time period, rather than each year
## e.g. if one hectad is recorded in multiple years within one time period, it would be included multiple times in the 10 lowest hectads
nmrsdata_wh_gs2 <- ddply(nmrsdata_wh_gs, .(Grid_square,Common_Name,easting,northing,Time_period,elevation1x1km,lat,lon), summarise,
                          PRESENCE = 1) 
## now order elevation by ascending value
## and take the top 10 values 
low_10_elev_wh_gs <- nmrsdata_wh_gs2 %>% group_by(Common_Name, Time_period) %>% arrange(elevation1x1km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_wh_gs <- summarySE(low_10_elev_wh_gs, measurevar="elevation1x1km", groupvars=c("Common_Name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_wh_gs, 
          elevation1x1km[Time_period == "TP1"] - elevation1x1km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.8782 - distribution is normal
qqnorm(d)
qqline(d) ## ok

## wilcoxon signed-rank test instead
low_elev_wh <- t.test(elevation1x1km ~ Time_period, data = low_elev_wh_gs, paired = TRUE)
low_elev_wh ## non-significant

## plot anyway
ggpaired(low_elev_wh_gs, x = "Time_period", y = "elevation1x1km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="t.test", paired = TRUE)

########################################
## Latitudinal shifts (trailing edge boundary)

## find 10 most southerly hectads & calculate mean latitude across years

## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_wh_gs <- nmrsdata_wh_gs2 %>% group_by(Common_Name, Time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_wh_gs <- summarySE(low_10_lat_wh_gs, measurevar="lat", groupvars=c("Common_Name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_lat_wh_gs, 
          lat[Time_period == "TP1"] - lat[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.07811 - normal distribution
qqnorm(d)
qqline(d) ## looks ok..

## wilcoxon signed-rank test instead
low_lat_wh <- t.test(lat ~ Time_period, data = low_lat_wh_gs, paired = TRUE)
low_lat_wh ## non-significant

## plot anyway
ggpaired(low_lat_wh_gs, x = "Time_period", y = "lat",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="t.test", paired = TRUE)
## very little change over time


########################################################
## heavily recorded grid squares
rm(list = ls())

nmrsdata_heavy_gs <- readRDS("../../Data/NMRS/NMRSdata_heavy_gs_98_02_12_16.rds")

## Filter 1: upland moth species only (either obligate or predominate use of uplands)
trait_data <- read.csv("../../Data/ecological_traits.csv", header=TRUE)

length(unique(nmrsdata_heavy_gs$Common_Name)) ## 759 species 

## subset species by upland levels 1 and 2 
north_species <- subset(trait_data, X7_montane_upland == 1 | X7_montane_upland == 2)
## remove nymphalids (butterflies)
north_species <- north_species[!(north_species$family=="Nymphalidae"),] ## 36 species
north_species <- north_species[,c(3,4,179)]
#write.csv(north_species, "../../Data/north_species2.csv", row.names=F)

## subset species in nmrs data
nmrsdata_heavy_gs <- merge(nmrsdata_heavy_gs, north_species, by.x="Preferred_Name", by.y="scientific_name")
length(unique(nmrsdata_heavy_gs$Common_Name)) ## 24 species

####### Filter 2: species must have enough coverage to calculate range shifts
####### Option 2: species need to occupy at least 20 hectads in each time period
nmrsdata_heavy_gs$Time_period <- ifelse(nmrsdata_heavy_gs$Year==2012 | nmrsdata_heavy_gs$Year==2013 | nmrsdata_heavy_gs$Year==2014
                                     | nmrsdata_heavy_gs$Year==2015 | nmrsdata_heavy_gs$Year==2016, "TP2", "TP1")
# get rid of the year information, keeping only which hectads each species was recorded in:
nmrsdata_heavy_gs2 <- ddply(nmrsdata_heavy_gs, .(Grid_square,Preferred_Name,Common_Name,easting,northing,Time_period), summarise,
                         PRESENCE = 1) 
# and now count the number of hectads per species for each year
tp_species_gs <- ddply(nmrsdata_heavy_gs2, .(Preferred_Name,Common_Name, Time_period), summarise,
                       Grid_squares = sum(PRESENCE))
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_gs <- tp_species_gs %>% group_by(Common_Name) %>% filter(all(Grid_squares>=20))
species_keep <- unique(tp_species_gs$Common_Name) ## 2 species 
## mostly removes rare/restricted species
nmrsdata_heavy_gs <-  nmrsdata_heavy_gs[which(nmrsdata_heavy_gs$Common_Name %in% species_keep),] ## 5 species
nmrsdata_heavy_gs <- unique(nmrsdata_heavy_gs) ## 289
length(unique(nmrsdata_heavy_gs$Grid_square)) ## 86 grid squares

## map
worldmap = map_data('world')
heavy_gs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_heavy_gs, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
heavy_gs
ggsave(wh_gs, file="../../Maps/Upland_moth_2spp_heavy_gs_1km.png")

### now calculate centroid elevation shift between time periods
## centroid = average elevation of all sampled occurrences during each time period
## and measure the low elevation shift - mean elevation of 10 lowest elevation hectads for each time period

## mean elevation for each species and each time period + SE + confidence intervals
library(Rmisc)
mean_elev_wh_gs <- summarySE(nmrsdata_wh_gs, measurevar="elevation1x1km", groupvars=c("Common_Name","Time_period"))

## check for normality 
# compute the difference
d <- with(mean_elev_wh_gs, 
          elevation1x1km[Time_period == "TP1"] - elevation1x1km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.3069 - distribution of the differences are not significantly different from normal distribution 
## i.e. they follow a normal distribution so t-test can be used
qqnorm(d)
qqline(d) ## looks fine
## northern dart, Rannoch Looper, Weaver's Wave
# mean_elev_rec_gs <- mean_elev_rec_gs[!(mean_elev_rec_gs$Common_Name=="Northern Dart"),]
# mean_elev_rec_gs <- mean_elev_rec_gs[!(mean_elev_rec_gs$Common_Name=="Rannoch Looper"),]
# mean_elev_rec_gs <- mean_elev_rec_gs[!(mean_elev_rec_gs$Common_Name=="Weaver's Wave"),]

mean_elev_wh <- t.test(elevation1x1km ~ Time_period, data = mean_elev_wh_gs, paired = TRUE)
mean_elev_wh ## non-significant

## plot anyway
library(ggpubr)
ggpaired(mean_elev_wh_gs, x = "Time_period", y = "elevation1x1km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Mean elevation")+
  stat_compare_means(method="t.test", paired = TRUE)
## non-significant downwards shift..

#####################################
## low elevational boundary shift
## this is the 10 lowest OCCUPIED hectads - can't have the same hectad repeated in this 10 just because it was recorded in multiple years

## first remove year info - this means we can only look at occupied hectads within each time period, rather than each year
## e.g. if one hectad is recorded in multiple years within one time period, it would be included multiple times in the 10 lowest hectads
nmrsdata_wh_gs2 <- ddply(nmrsdata_wh_gs, .(Grid_square,Common_Name,easting,northing,Time_period,elevation1x1km,lat,lon), summarise,
                         PRESENCE = 1) 
## now order elevation by ascending value
## and take the top 10 values 
low_10_elev_wh_gs <- nmrsdata_wh_gs2 %>% group_by(Common_Name, Time_period) %>% arrange(elevation1x1km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_wh_gs <- summarySE(low_10_elev_wh_gs, measurevar="elevation1x1km", groupvars=c("Common_Name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_wh_gs, 
          elevation1x1km[Time_period == "TP1"] - elevation1x1km[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.8782 - distribution is normal
qqnorm(d)
qqline(d) ## ok

## wilcoxon signed-rank test instead
low_elev_wh <- t.test(elevation1x1km ~ Time_period, data = low_elev_wh_gs, paired = TRUE)
low_elev_wh ## non-significant

## plot anyway
ggpaired(low_elev_wh_gs, x = "Time_period", y = "elevation1x1km",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="t.test", paired = TRUE)

########################################
## Latitudinal shifts (trailing edge boundary)

## find 10 most southerly hectads & calculate mean latitude across years

## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_wh_gs <- nmrsdata_wh_gs2 %>% group_by(Common_Name, Time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_wh_gs <- summarySE(low_10_lat_wh_gs, measurevar="lat", groupvars=c("Common_Name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_lat_wh_gs, 
          lat[Time_period == "TP1"] - lat[Time_period == "TP2"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = 0.07811 - normal distribution
qqnorm(d)
qqline(d) ## looks ok..

## wilcoxon signed-rank test instead
low_lat_wh <- t.test(lat ~ Time_period, data = low_lat_wh_gs, paired = TRUE)
low_lat_wh ## non-significant

## plot anyway
ggpaired(low_lat_wh_gs, x = "Time_period", y = "lat",
         color = "Time_period", line.color = "gray", line.size = 0.4,
         palette = "jco", id="Common_Name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="t.test", paired = TRUE)
## very little change over time























