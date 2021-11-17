##########################
#### user: Lisbeth Hordley
#### date: October 2021
#### info: Filtering NMRS data with elevation & temperature to cool-adapted moths

rm(list = ls())

## load packages
library(dplyr)
library(ggplot2)

## NMRS temperature data for each TP
nmrsdata_early <- read.csv("Data/NMRS/NMRS_hectad_elevation_climate_TP1.csv", header=TRUE) ## NMRS data for TP1 with elevation and climate data 
nmrsdata_late <- read.csv("Data/NMRS/NMRS_hectad_elevation_climate_TP2.csv", header=TRUE) ## NMRS data for TP2 with elevation and climate data 
## migrant hectads that need excluded
migrant_hectads_early <- read.csv("Data/NMRS/migrant_hectad_exclusion_TP1.csv", header=TRUE) ## species x hectad combinations to remove due to immigrant populations
migrant_exclusion_late <- read.csv("Data/NMRS/migrant_hectad_exclusion_TP2.csv", header=TRUE) ## both provided by Richard on Teams (20/10/2021)
## upland species list to filter
cool_species <- read.csv("Data/NMRS/NMRS_temp_median_percentile_selection.csv", header=TRUE) ## provided by Richard on Teams (20/10/2021)

### first exclude migrant hectads from each nmrs dataset
migrant_hectads_early[is.na(migrant_hectads_early)] <- 0
nmrsdata_early <- merge(nmrsdata_early, migrant_hectads_early, by=c("Common_Name", "Hectad"), all.x=TRUE)
nmrsdata_early[is.na(nmrsdata_early)] <- 0
nmrsdata_early <- nmrsdata_early[nmrsdata_early$Exclude !=1, ] 

migrant_exclusion_late[is.na(migrant_exclusion_late)] <- 0
nmrsdata_late <- merge(nmrsdata_late, migrant_exclusion_late, by=c("Common_Name", "Hectad"), all.x=TRUE)
nmrsdata_late[is.na(nmrsdata_late)] <- 0
nmrsdata_late <- nmrsdata_late[nmrsdata_late$Exclude !=1, ] 

## filter by cool adapted species only 
cool_species[is.na(cool_species)] <- 0 # replace NAs with 0 (0 = include the species as northern/upland)
cool_species_final <- cool_species[cool_species$exclude !=1, ] ## 73 species
cool_sp <- cool_species_final[,1]
nmrsdata_early <- nmrsdata_early[which(nmrsdata_early$Common_Name %in% cool_sp), ] # 5361 rows
nmrsdata_late <- nmrsdata_late[which(nmrsdata_late$Common_Name %in% cool_sp), ] # 11013 rows
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
nmrsdata_final <- rbind(nmrsdata_early, nmrsdata_late) ## 16373 rows
length(unique(nmrsdata_final$Hectad)) ## 1976 unique hectads
## save file
saveRDS(nmrsdata_final, file="Data/NMRS/NMRS_cool_moths_final.rds")

