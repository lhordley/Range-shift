##########################
#### user: Lisbeth Hordley
#### date: October 2021
#### info: Filtering NMRS data with elevation & temperature to cool-adapted moths

rm(list = ls())

## load packages
library(dplyr)
library(ggplot2)

## NMRS data with climate and elevation
nmrsdata <- readRDS("Data/NMRS/NMRS_hectad_elevation_climate.rds") ## 667,521 rows
## migrant hectads that need excluded
migrant_hectads <- read.csv("Data/NMRS/NMRS_migrant_hectads_exclude.csv", header=TRUE) ## species x hectad combinations to remove due to immigrant populations
## upland species list to filter
cool_species <- read.csv("Data/NMRS/NMRS_new_temp_median_percentile_selection.csv", header=TRUE) ## provided by Richard on Teams (16/02/2022)

## Chimney Sweeper name is incorrect in NMRS - change now
nmrsdata$Common_name[nmrsdata$Common_name == "Chimney Sweep"] <- "Chimney Sweeper"

### first exclude migrant hectads from each nmrs dataset
migrant_hectads[is.na(migrant_hectads)] <- 0
nmrsdata <- merge(nmrsdata, migrant_hectads, by=c("Common_name", "Hectad", "Time_period"), all.x=TRUE)
nmrsdata[is.na(nmrsdata)] <- 0
nmrsdata <- nmrsdata[nmrsdata$Exclude !=1, ] ## 667,288 rows

## filter by cool adapted species only 
cool_species[is.na(cool_species)] <- 0 # replace NAs with 0 (0 = include the species as northern/upland)
cool_species_final <- cool_species[cool_species$Exclude !=1, ] ## 77 species
cool_sp <- cool_species_final[,1]
nmrsdata <- nmrsdata[which(nmrsdata$Common_name %in% cool_sp), ] # 19230 rows
length(unique(nmrsdata$Common_name)) # 77 species 

sci_name <- unique(nmrsdata[,c("Common_name","Scientific_name")])
write.csv(sci_name, file="sci_names.csv", row.names=FALSE)

## check that each species is found in both time periods
table(nmrsdata$Common_name, nmrsdata$Time_period)
## all 77 species are found in both time periods

## save file
saveRDS(nmrsdata, file="Data/NMRS/NMRS_cool_moths_final.rds")

## Make maps of migrant hectads
lat_lon <- unique(nmrsdata[,c("lat","lon","Hectad")])
migrant_hectads <- merge(migrant_hectads, lat_lon, by="Hectad", all.x=TRUE)
migrant_hectads$Exclude <- as.factor(migrant_hectads$Exclude)
worldmap = map_data('world')
migrant_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = migrant_hectads, 
             aes(x = as.numeric(lon), 
                 y = as.numeric(lat), colour=Exclude)) + 
  scale_color_manual(values = c("black","red")) +
  theme_void() +
  theme(title = element_text(size = 12)) +
  theme(legend.position = "none") +
  facet_wrap(~ Common_name, ncol=4, nrow=2) +
  theme(strip.text.x = element_text(size = 12))
migrant_hecs

ggsave(migrant_hecs, file="Outputs/Maps/FigureS1_immigrant_hecs.png", height=10, width=15)


##############################################################################################

## Look at distribution trends for cool moths

trends <- read.csv("Data/Species trends for Moth Atlas (final).csv", header=TRUE)
cool_species <- read.csv("Data/NMRS/NMRS_new_temp_median_percentile_selection.csv", header=TRUE) ## provided by Richard on Teams (16/02/2022)

# Remove capital B in trends so it matches NMRS names
trends$Vernacular[trends$Vernacular == "Manchester Treble-Bar"] <- "Manchester Treble-bar"
# Remove Lunar Thorn
trends <- trends[!trends$Vernacular=="Lunar Thorn",]
cool_species <- cool_species[!cool_species$Common_name=="Lunar Thorn",]

cool_species[is.na(cool_species)] <- 0 # replace NAs with 0 (0 = include the species as northern/upland)
cool_species_final <- cool_species[cool_species$Exclude !=1, ] ## 76 species
cool_sp <- cool_species_final[,1]


# Long term trend stats: 1970-2016
trends_cool <- trends[which(trends$Vernacular %in% cool_sp), ] # 76 species
trends_cool <- trends_cool %>% mutate_all(na_if,"") # change blanks to NA
trends_cool <- trends_cool[,c("Vernacular", "TREND.DLT", "SIG.DLT", "DATE.DLT")]

colSums(is.na(trends_cool)) # 44 species with NAs - could not calculate a trend
trends_cool <- na.omit(trends_cool) # 32 species with trends

str(trends_cool)

# calculate 4 stats:
# 1. Average % of distribution change for all species
trends_cool$TREND.DLT <- gsub("%", "", trends_cool$TREND.DLT)
trends_cool$TREND.DLT <- as.numeric(trends_cool$TREND.DLT)
mean(trends_cool$TREND.DLT) # -33.25% decline across all 32 species

# 2. Average % of distribution change for SIGNIFICANTLY changing species
mean(trends_cool[trends_cool$SIG.DLT == 'Y', 'TREND.DLT']) # -38.96% significant decline across all 32 species

# 3. % of species which are declining
nrow(trends_cool[trends_cool$TREND.DLT<0,]) # 27 species = 84.38%

# 4. % of species which are SIGNIFICANTLY declining
nrow(trends_cool[trends_cool$TREND.DLT<0 & trends_cool$SIG.DLT=="Y",]) # 23 species = 71.88%

# make graph to illustrate?
hist(trends_cool$TREND.DLT) # could create this in ggplot

## Short term trend stats: 2000-2016
trends_cool2 <- trends[which(trends$Vernacular %in% cool_sp), ] # 76 species
trends_cool2 <- trends_cool2 %>% mutate_all(na_if,"") # change blanks to NA
trends_cool2 <- trends_cool2[,c("Vernacular", "TREND.DST", "SIG.DST", "DATE.DST")]
colSums(is.na(trends_cool2)) # 38 species with NAs - could not calculate a trend
trends_cool2 <- na.omit(trends_cool2) # 38 species with trends

# calculate 4 stats:
# 1. Average % of distribution change for all species
trends_cool2$TREND.DST <- gsub("%", "", trends_cool2$TREND.DST)
trends_cool2$TREND.DST <- as.numeric(trends_cool2$TREND.DST)
mean(trends_cool2$TREND.DST) # -9.53% decline across all 38 species

# 2. Average % of distribution change for SIGNIFICANTLY declining species
mean(trends_cool2[trends_cool2$SIG.DST == 'Y', 'TREND.DST']) # -24% significant decline across all 38 species

# 3. % of species which are TREND.DST
nrow(trends_cool2[trends_cool2$TREND.DST<0,]) # 27 species = 71.05%

# 4. % of species which are SIGNIFICANTLY declining
nrow(trends_cool2[trends_cool2$TREND.DST<0 & trends_cool2$SIG.DST=="Y",]) # 5 species = 13.16%

