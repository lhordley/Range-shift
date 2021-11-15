##########################
#### user: Lisbeth Hordley
#### date: October 2021
#### info: Filtering NMRS data with elevation & temperature to cool-adapted moths

rm(list = ls())

## load packages
library(dplyr)
library(ggplot2)

## NMRS temperature data for each TP
nmrsdata_temp_early <- readRDS("Data/NMRS/NMRS_hectad_temperature_TP1.rds") ## NMRS data for TP1 with temperature 
nmrsdata_temp_late <- readRDS("Data/NMRS/NMRS_hectad_temperature_TP2.rds") ## NMRS data for TP2 with temperature 
## migrant hectads that need excluded
migrant_hectads_early <- read.csv("Data/NMRS/migrant_hectad_exclusion_TP1.csv", header=TRUE) ## species x hectad combinations to remove due to immigrant populations
migrant_exclusion_late <- read.csv("Data/NMRS/migrant_hectad_exclusion_TP2.csv", header=TRUE) ## both provided by Richard on Teams (20/10/2021)
## upland species list to filter
cool_species <- read.csv("Data/NMRS/NMRS_temp_median_percentile_selection.csv", header=TRUE) ## provided by Richard on Teams (20/10/2021)

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

## filter by cool adapted species only 
cool_species[is.na(cool_species)] <- 0 # replace NAs with 0 (0 = include the species as northern/upland)
cool_species_final <- cool_species[cool_species$exclude !=1, ] ## 73 species
cool_sp <- cool_species_final[,1]
nmrsdata_temp_early <- nmrsdata_temp_early[which(nmrsdata_temp_early$Common_Name %in% cool_sp), ] # 5390 rows
nmrsdata_temp_late <- nmrsdata_temp_late[which(nmrsdata_temp_late$Common_Name %in% cool_sp), ] # 11079 rows
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
nmrsdata_temp_final <- rbind(nmrsdata_temp_early, nmrsdata_temp_late) ## 16373 rows
length(unique(nmrsdata_temp_final$Hectad)) ## 1976
## remove centred coordinates and exclude columns
nmrsdata_temp_final[,c("lat_centre","lon_centre","easting_centre","northing_centre","Exclude"):=NULL] 
## save file
saveRDS(nmrsdata_temp_final, file="Data/NMRS/NMRS_cool_moths_final.rds")

