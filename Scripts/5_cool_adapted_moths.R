##########################
#### user: Lisbeth Hordley
#### date: October 2021
#### info: Filtering NMRS data with elevation & temperature to cool-adapted moths

rm(list = ls())

## load packages
library(dplyr)
library(ggplot2)

########################################################################################################
## STEP 1: PLOT SPECIES TEMPERATURE RANGES (MEDIAN & 75TH PERCENTILE) BETWEEN 1975-1991
## These are given to Richard to determine which cool-adapted species to select
nmrsdata_climate <- readRDS("Data/NMRS/NMRS_hectad_elevation_climate.rds") ## NMRS data with temperature 
migrant_hectads <- read.csv("Data/NMRS/NMRS_migrant_hectads_exclude.csv", header=TRUE) ## species x hectad combinations to remove due to immigrant populations
## provided by Richard in Teams 20/10/21

## Remove TP1 data - only looking at temperature range between 1975 - 1991
nmrsdata_climate <- nmrsdata_climate[nmrsdata_climate$Time_period=="1975-1991",] # 254689
migrant_hectads <- migrant_hectads[migrant_hectads$Time_period=="1975-1991",]
## merge migrant exclusion data with nmrs_temp data
## remove hectad/species combinations where exclude == 1
migrant_hectads[is.na(migrant_hectads)] <- 0
nmrsdata_climate <- merge(nmrsdata_climate, migrant_hectads, by=c("Common_name", "Hectad", "Time_period"), all=TRUE)
nmrsdata_climate[is.na(nmrsdata_climate)] <- 0
nmrsdata_climate <- nmrsdata_climate[nmrsdata_climate$Exclude !=1, ] ## 254544 rows

detach(package:plyr)
library(dplyr)
library(ggplot2)
summary_nmrs_temp <- nmrsdata_climate %>% 
  group_by(Common_name) %>% 
  summarise(lower = min(temperature), upper = max(temperature), p = median(temperature)) ## 788 species
summary_nmrs_temp$lower <- as.numeric(summary_nmrs_temp$lower)
summary_nmrs_temp$upper <- as.numeric(summary_nmrs_temp$upper)
summary_nmrs_temp$p <- as.numeric(summary_nmrs_temp$p)

summary_nmrs_temp <- summary_nmrs_temp[order(summary_nmrs_temp$p, decreasing = TRUE),]  
library(Hmisc)
summary_nmrs_temp$groups<-as.numeric(cut2(summary_nmrs_temp$p, g=12))
## these groups are used to put groups into nmrs_temp to plot them more easily
sp_groups <- summary_nmrs_temp[,c(1,5)]

library(Hmisc)
## use groups from above to put each species into groups to make plots
nmrsdata_climate <- merge(nmrsdata_climate, sp_groups, by="Common_name", all=TRUE)
group <- unique(nmrsdata_climate$groups) ## 12 groups

for(i in group) {
  print(i)
  # Printing ggplot within for-loop
  temp_plot <- ggplot(data = nmrsdata_climate[nmrsdata_climate$groups==i,], mapping = aes(x = reorder(Common_name,-temperature, FUN=median), y = temperature)) +
    geom_boxplot() +
    #geom_hline(yintercept=9, linetype="dotted") +
    coord_flip() +
    labs(x="Common name", y="Median temperature") +
    scale_y_continuous(breaks = seq(3, 11, by = 1)) +
    theme_light()
  
  ggsave(temp_plot, file=paste0("Data/UKCP_tas_annual/Early_TP_75_91/NMRS_new_temp_range_box_", i,".png"), width = 20, height = 25, units = "cm")
  Sys.sleep(2)
}

## calculate median and 75th percentile for each species for Richard
temp_ranges <- nmrsdata_climate %>%  group_by(Common_name) %>% 
  summarise(median_temperature = median(temperature), percentile_75th = quantile(temperature, probs = 0.75))
## order by median ascending
temp_ranges <- arrange(temp_ranges, median_temperature)
## save file
write.csv(temp_ranges, file="Data/NMRS_new_temp_median_percentile.csv", row.names=FALSE)
## Boxplots and percentile list given to Richard to select cool adapted species

########################################################################################################
## STEP 2: FILTER NMRS DATA TO SELECTION OF COOL-ADAPTED SPECIES DEFINED BY RICHARD USING 75TH TEMPERATURE PERCETILE

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

