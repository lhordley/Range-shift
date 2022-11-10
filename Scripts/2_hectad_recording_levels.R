##########################
#### user: Lisbeth Hordley
#### date: July 2021
#### info: Hectad recording levels NMRS data

rm(list = ls())
options(scipen=999)

nmrsdata <- readRDS("Data/NMRS/NMRS_hectad_cleaned.rds") # 667,521 rows
## easting/northing, lat/lon and elevation are at 10km resolution

# count the number of recorded hectads and the number of species*hectad records in each year
library(plyr)
nmrsdata$PRESENCE <- 1
nmrs_rec_hectad_years <- ddply(nmrsdata, .(Hectad, Time_period, easting, northing), numcolwise(mean)) ## removes species data, just hectads recorded for each year

# Calculate recording effort
# This is the number of records for each hectad summed across T1 and T2
# To be used in extinction models 
recording_effort <- nmrsdata %>% group_by(Hectad) %>% summarise(n_records=n())
write.csv(recording_effort, file="Data/Recording_effort_all_spp.csv", row.names=FALSE)

# # let's actually count the number of recorded hectads and the number of species*hectad records in each year
# years.extended <- 1970:2016
# 
# recording <- data.frame(YEAR = numeric(),
#                         HECTADS = numeric(),
#                         RECORDS = numeric())
# 
# for (n in years.extended){
#   year.hec <- nmrs_rec_hectad_years[which(nmrs_rec_hectad_years$Year==n), ]
#   HECTADS <- nrow(year.hec)
# 
#   year.rec <- nmrsdata[which(nmrsdata$Year==n), ]
#   RECORDS <- nrow(year.rec)
# 
#   YEAR <- n
# 
#   out <- data.frame(cbind(YEAR,HECTADS,RECORDS))
#   recording <- rbind(recording,out)
# 
# 
# }
# 
# summary(recording)
# 
# 
# # now plot the number of records over time
# records_time <- ggplot(recording, aes(x=YEAR, y=RECORDS)) +
#   geom_point(size=2) +
#   labs(x="Year", y="Number of records") +
#   scale_y_continuous(breaks = seq(60000, 1400000, by=400000)) +
#   theme_classic()
# records_time
# 
# plot(RECORDS ~ YEAR, data = recording)

## time periods of interest
tp1 <- 1975:1991
tp2 <- 2012:2016

# first follow Hickling et al directly by calculating which grid squares were recorded at least once in each interval
# first select out the intervals
nmrs_rec_hectad_years_tp1 <- nmrs_rec_hectad_years[which(nmrs_rec_hectad_years$Time_period=="1975-1991"), ] # 2239
nmrs_rec_hectad_years_tp2 <- nmrs_rec_hectad_years[which(nmrs_rec_hectad_years$Time_period=="2012-2016"), ] # 2647
## 2239 in tp1
## 2647 in tp2

# plot these on a map
worldmap = map_data('world')
ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrs_rec_hectad_years_tp2,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat)), size=1, colour="blue") +
  geom_point(data =nmrs_rec_hectad_years_tp1,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat)), size=1, colour="red") +
  theme(title = element_text(size = 12))
## historical TP doesn't have great coverage in Scotland where we need it
## going back to 1970 doesn't make much difference - only adds a handful of grid squares in Scotland
## tried extending to 1995 too - also doesn't make a huge difference

## now we want to recombine these datasets in a way that ditches hectads not present in both TPs
nmrs_rec_hectad_tp1_vec <- as.character(nmrs_rec_hectad_years_tp1$Hectad)
nmrs_rec_hectad_tp2_vec <- as.character(nmrs_rec_hectad_years_tp2$Hectad)

nmrs_rec_hectad_full <- append(nmrs_rec_hectad_tp1_vec,nmrs_rec_hectad_tp2_vec)

nmrs_rec_hectad_full <- data.frame(cbind(nmrs_rec_hectad_full,1)) ## list of all grid squares in both TPs
colnames(nmrs_rec_hectad_full) <- c("Hectad","Intervals")
nmrs_rec_hectad_full$Hectad <- as.factor(as.character(nmrs_rec_hectad_full$Hectad)) ## 4886
nmrs_rec_hectad_full$Intervals <- as.numeric(as.character(nmrs_rec_hectad_full$Intervals))


nmrs_rec_hectad_good <- ddply(nmrs_rec_hectad_full, .(Hectad), numcolwise(sum)) ## hectads either have 1 (only recorded in one TP), or 2 (recorded in both TPs)

nmrs_rec_hectad_good <- nmrs_rec_hectad_good[which(nmrs_rec_hectad_good$Intervals > 1), ] ## only keep hectads recorded in both (i.e. intervals > 1)
## 2191 good hectads

# drop unused levels
nmrs_rec_hectad_good <- droplevels(nmrs_rec_hectad_good)

# generate the vector
nmrs_good_hectads <- levels(nmrs_rec_hectad_good$Hectad)

# pick out the good records
length(unique(nmrsdata$Hectad)) ## 2695 total hectads
nmrsdata_good_hectads <- nmrsdata[which(nmrsdata$Hectad %in% nmrs_good_hectads), ]
length(unique(nmrsdata_good_hectads$Hectad)) ## 2191
# we lose ~20% of hectads

nmrs_good_hectad2 <- ddply(nmrsdata_good_hectads, .(Hectad, Time_period), numcolwise(mean))

plot <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrs_good_hectad2, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat)), size=1, colour="purple") +
  theme(title = element_text(size = 12))
plot ## these are recorded hectads (2191 total)


## we can also further restrict to the hectads that have records of a minimum percentage of the regional species richness
# this is going to take quite a complex loop to construct:
# for every hectad in the NMRS_GB_good dataframe, we need to calculate the observed species richness
# then we need to calculate pairwise distances to all other recorded hectads, pick out the nearest 100, 
# and calculate the regional species richness of these 100
# and finally turn this into a percentage of regional SR that is recorded in the hectad

# start with the list of recorded hectads, re-attach the eastings and northings
nmrs_good <- ddply(nmrs_rec_hectad_years, .(Hectad,easting,northing), summarise,
                   COUNT = sum(PRESENCE))

# select out only the good hectads
nmrs_good <- nmrs_good[which(nmrs_good$Hectad %in% nmrs_good_hectads), ] # 2191
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
hec_records <- data.frame(Time_period = factor(),
                          Hectad = factor(),
                          Easting = numeric(),
                          Northing = numeric(),
                          hec_SR = numeric(),
                          reg_SR = numeric(),
                          hec_perc = numeric())
i <- 1

for (x in nmrs_good_hectads){ # take each good hectad
  if (round(i, -1) == i){
    print(i)
  }
  hec <- nmrs_good[which(nmrs_good$Hectad == x), ]    # pick out details of focal hectad
  candidates <- nmrs_good[which(nmrs_good$Hectad != x), ]   # pick out details of all others (2190 hectads)
  
  hec_east <- hec[[1,2]]     # easting of target
  hec_north <- hec[[1,3]]    # northing of target
  
  candidates$east_diff <- candidates$easting - hec_east             # longitudinal difference
  candidates$north_diff <- candidates$northing - hec_north          # latitudinal difference
  
  candidates$distance <- sqrt((candidates$east_diff^2) + (candidates$north_diff^2))    # absolute difference
  
  candidates <- candidates[order(candidates$distance),] # sort by distance ascending
  
  closest <- candidates[1:100,] # select out closest 100
  
  # now we want to calculate species richness from the hectad and from the region in both time periods
  for (n in unique(nmrsdata_good_hectads$Time_period)){
    year_recs <- nmrsdata_good_hectads[nmrsdata_good_hectads$Time_period==n, ] ## nmrs good hectad data for tp of interest
    
    hec_recs <- year_recs[which(year_recs$Hectad == x), ] # pull out hectad records
    hec_SR <- length(unique(hec_recs$Common_name))       # calculate species richness
    
    reg_recs <- year_recs[which(year_recs$Hectad %in% closest$Hectad), ]  # pull out region records
    reg_recs <- rbind(reg_recs,hec_recs)    # add in hectad records (they're part of the regional richness too!)
    reg_SR <- length(unique(reg_recs$Common_name))
    
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
lat_lon <- nmrsdata_good_hectads %>% distinct(Hectad, lat, lon, .keep_all = FALSE)

hec_records <- merge(hec_records, lat_lon, by.x="HECTAD", by.y="Hectad")

# this data takes a while to generate so let's back it up
write.csv(hec_records, "Data/NMRS/Hectad_recording_levels_1975_1991_2012_2016.csv", row.names = F)

# ## Analysis 1: hectads with at least one species recorded in BOTH time periods: 1782 hectads
# ## file = nmrsdata_good_hectads (needs subsetting by years)
# years <- c(1975:1991, 2012:2016)
# nmrsdata_rec_hecs <- nmrsdata_good_hectads[which(nmrsdata_good_hectads$Year %in% years), ]
# length(unique(nmrsdata_rec_hecs$Hectad)) ## 1782 hectads
# length(unique(nmrsdata_rec_hecs$Year)) ## 22 years
# ## save file
# saveRDS(nmrsdata_rec_hecs, "../../Data/NMRS/NMRSdata_rec_hecs_75_91_12_16.rds")
# 
# ## Analysis 2a: hectads classed as well-recorded in BOTH time periods: 57 hectads (not doable)
# well_hecs <- hec_records %>%
#   filter(RECORDING.LEVEL == "Well recorded") %>%
#   group_by(HECTAD) %>%
#   dplyr::summarise(n_row=n())
# well_hecs <- well_hecs[which(well_hecs$n_row > 1), ] 
# well_hecs$n_row <- NULL
# 
# ## Analysis 2b: hectads classed as heavily recorded in BOTH time periods: 412 hectads
# heavy_hecs <- hec_records %>%
#   filter(RECORDING.LEVEL == "Heavily recorded") %>%
#   group_by(HECTAD) %>%
#   dplyr::summarise(n_row=n())
# heavy_hecs <- heavy_hecs[which(heavy_hecs$n_row > 1), ] 
# heavy_hecs$n_row <- NULL
# ## merge into nmrs data and filter by years so only have years within time periods
# nmrsdata_heavy_hecs <- nmrsdata[which(nmrsdata$Hectad %in% heavy_hecs$HECTAD)]
# length(unique(nmrsdata_heavy_hecs$Hectad)) ## 412 hectads
# nmrsdata_heavy_hecs <- nmrsdata_heavy_hecs[which(nmrsdata_heavy_hecs$Year %in% years)]
# length(unique(nmrsdata_heavy_hecs$Year)) ## 22 years
# ## save file
# saveRDS(nmrsdata_heavy_hecs, "../../Data/NMRS/NMRSdata_heavy_hecs_75_91_12_16.rds")
# 
# ## Analysis 2c: hectads classed as well OR heavily recorded in BOTH time periods: 747 hectads
# well_heavy_hecs <- hec_records %>%
#   filter(RECORDING.LEVEL == "Heavily recorded" | RECORDING.LEVEL=="Well recorded") %>%
#   group_by(HECTAD) %>%
#   dplyr::summarise(n_row=n())
# well_heavy_hecs <- well_heavy_hecs[which(well_heavy_hecs$n_row > 1), ] 
# well_heavy_hecs$n_row <- NULL
# ## merge into nmrs data and filter by years so only have years within time periods
# nmrsdata_well_heavy_hecs <- nmrsdata[which(nmrsdata$Hectad %in% well_heavy_hecs$HECTAD)]
# length(unique(nmrsdata_well_heavy_hecs$Hectad)) ## 747 hectads
# nmrsdata_well_heavy_hecs <- nmrsdata_well_heavy_hecs[which(nmrsdata_well_heavy_hecs$Year %in% years)]
# length(unique(nmrsdata_well_heavy_hecs$Year)) ## 22 years
# ## save file
# saveRDS(nmrsdata_well_heavy_hecs, "../../Data/NMRS/NMRSdata_well_heavy_hecs_75_91_12_16.rds")
