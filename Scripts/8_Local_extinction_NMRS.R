##########################
#### user: Lisbeth Hordley
#### date: December 2021
#### info: Drivers of local extinction

rm(list = ls())

library(ggplot2)
library(ggeffects)
library(DHARMa)
library(lme4)
library(MuMIn)
library(stringr)
library(blmeco)
library(MuMIn)
library(dplyr)
library(ggpubr)
library(data.table)
library(tidyr)
options(scipen=999)

########################## Step 1: Get data into correct format for modelling ########################## 

#####################################
## well & heavily recorded hectads ##
#####################################

# # First calculate recording effort - read in data on all NMRS records for years of interest (1975-1991 and 2012-2016)
# tot_records <- read.csv("Data (all)/lhordley_export.csv", header=TRUE)
# # Calculate recording effort to put in future models
# # This is calculated as the total number of records in T1 and T2 (summed together)
# recording_effort <- tot_records %>% group_by(Hectad) %>% summarise(n_recs=sum(nRecords))
# write.csv(recording_effort, file="Data/Recording_effort.csv", row.names=FALSE)
# 
# # Number of records in T1 and T2 for response to reviewers (based on well recorded hectads only)
# well_recs <- unique(nmrsdata_well$Hectad)
# tot_records_well <- tot_records[which(tot_records$Hectad %in% well_recs), ]
# recording_effort2 <- tot_records %>% group_by(Year) %>% summarise(n_recs=sum(nRecords))
# recording_effort2$TP <- ifelse(recording_effort2$Year<=1991, "TP1", "TP2")
# recording_effort3 <- recording_effort2 %>% group_by(TP) %>% summarise(mean=mean(n_recs))

## Read in data
# NMRS data for cool-adapted moths with elevation and temperature data at a 10km scale
nmrsdata <- readRDS("Data/NMRS_cool_moths_final.rds") ## NMRS data for all hectads and all years with elevation
## hectad recording levels
hec_records <- read.csv("Data/Hectad_recording_levels_1975_1991_2012_2016.csv", header=TRUE)
# recording effort
recording_effort <- read.csv("Data/Recording_effort.csv", header=TRUE)

nmrsdata <- nmrsdata[!nmrsdata$Common_name=="Lunar Thorn",]

## Step 1a. Filter NMRS data to well or heavily recorded hectads only
well_heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded" | RECORDING.LEVEL=="Well recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
well_heavy_hecs <- well_heavy_hecs[which(well_heavy_hecs$n_row > 1), ]
well_heavy_hecs$n_row <- NULL ## 1084 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_well <- nmrsdata[which(nmrsdata$Hectad %in% well_heavy_hecs$HECTAD), ]
length(unique(nmrsdata_well$Hectad)) ## 962 hectads
length(unique(nmrsdata_well$Common_name)) ## 76 species

## Step 1b. Define colonised, persisted, and extirpated hectads
nmrsdata_well$Recorded <- 1
nmrsdata_well_expand <- nmrsdata_well %>% tidyr::expand(Common_name, Time_period, Hectad)
recorded <- nmrsdata_well[,c("Common_name", "Time_period", "Hectad", "Recorded")]
nmrsdata_well_expand <- merge(nmrsdata_well_expand, recorded, by=c("Common_name", "Time_period", "Hectad"), all.x=TRUE)
nmrsdata_well_expand[is.na(nmrsdata_well_expand)] <- 0 ## change NAs to zero == species was NOT recorded at this hectad in this time period
## change to long format
nmrsdata_well_expand2 <- nmrsdata_well_expand %>%
  spread(Time_period, Recorded) ## each species has 1424 rows = the number of recorded hectads
## first remove rows where TP1 AND TP2 == 0 (this a hectad where a species was never recorded)
nmrsdata_well_expand2<-nmrsdata_well_expand2[!(nmrsdata_well_expand2$`1975-1991`==0 & nmrsdata_well_expand2$`2012-2016`==0),]
nmrsdata_well_expand2$Hectad_category <- case_when(
  nmrsdata_well_expand2$`1975-1991`==0 & nmrsdata_well_expand2$`2012-2016`==1 ~ "Colonisation",
  nmrsdata_well_expand2$`1975-1991`==1 & nmrsdata_well_expand2$`2012-2016`==0 ~ "Extirpation",
  TRUE ~ "Persistence"
)

## Step 1c. Plot the proportion of extinctions at each hectad on a map
## add lat/lon data back in
lat_lon <- nmrsdata_well[,c("Hectad", "lat", "lon","elevation10x10km","elevation10x10km_SD")]
lat_lon <- lat_lon %>% distinct(Hectad, .keep_all = TRUE) ## different way of unique
nmrsdata_well_expand2 <- merge(nmrsdata_well_expand2, lat_lon, by="Hectad", all.x=TRUE)

## find proportion of species at each extirpated hectad - where are extirpations occurring? 
extirpated_hecs1 <- nmrsdata_well_expand2 %>% group_by(Hectad, lat, lon) %>% filter(`1975-1991`==1) %>%
  dplyr::summarise(tot_sp = n()) # number of species recorded in each hectad in TP1
extirpated_hecs2 <- nmrsdata_well_expand2 %>% group_by(Hectad, lat, lon) %>% filter(sum(`1975-1991`)==0) %>%
  dplyr::summarise(tot_sp = 0)
extirpated_hecs3 <- nmrsdata_well_expand2 %>% group_by(Hectad, lat, lon) %>%
  filter(Hectad_category=="Extirpation") %>%
  dplyr::summarise(extir_sp = n()) # number of species which have gone extinct from each hectad
extirpated_hecs <- rbind(extirpated_hecs1, extirpated_hecs2)
extirpated_hecs <- merge(extirpated_hecs, extirpated_hecs3, by=c("Hectad", "lat", "lon"), all.x=TRUE)
extirpated_hecs[is.na(extirpated_hecs)] <- 0 ## change NAs to zero - where there are no extirpations at a hectad
extirpated_hecs$extir_prop <- extirpated_hecs$extir_sp/extirpated_hecs$tot_sp ## 962
extirpated_hecs$extir_prop[is.nan(extirpated_hecs$extir_prop)]<-0 # convert NaNs to zero (hectads where there are only colonising species - the 93 hectads in extirpated_hecs2)

## map
worldmap = map_data('world')
well_extir_hecs <- ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  theme_void() +
  geom_point(data =extirpated_hecs,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat), colour=extir_prop), size=1) +
  scale_color_viridis_c(name="Proportion of species") + 
  theme(title = element_text(size = 12)) +
  scale_fill_continuous(guide = guide_legend()) +
  guides(colour = guide_colourbar(title.vjust = 0.9)) +
  theme(legend.position="bottom", legend.key.height = unit(0.5, "cm"), legend.key.width = unit(2, "cm"))
well_extir_hecs

## Step 1d. Add in recording effort and climate data and calculate difference between T1 and T2 temperature and precipitation
# add recording effort in
nmrsdata_well_expand2 <- merge(nmrsdata_well_expand2, recording_effort, by="Hectad", all.x=TRUE)
## add in climate data to each data frame
nmrs_climate <- readRDS("Data/NMRS_hectad_elevation_climate.rds")
climate_hecs <- unique(nmrs_climate[,c(1,2,16,17)]) # 2719 rows
# change from long to wide
climate_wide <- reshape(climate_hecs, idvar = "Hectad", timevar = "Time_period", direction = "wide")
well_recs <- unique(nmrsdata_well$Hectad)
climate_wide2 <- climate_wide[which(climate_wide$Hectad %in% well_recs), ]
colnames(climate_wide2)[2:5] <- c("temperature.TP1", "total_precip.TP1", "temperature.TP2", "total_precip.TP2")
climate_wide2$temp_diff <- climate_wide2$temperature.TP2 - climate_wide2$temperature.TP1
climate_wide2$precip_diff <- climate_wide2$total_precip.TP2 - climate_wide2$total_precip.TP1
## merge in with nmrs data
nmrsdata_well_expand2 <- merge(nmrsdata_well_expand2, climate_wide2, by="Hectad", all.x=TRUE)

## Step 1e. Remove colonisations - only interested in extinctions and persistence 
nmrsdata_well_expand2<-nmrsdata_well_expand2[!(nmrsdata_well_expand2$Hectad_category=="Colonisation"),] ## remove colonisations
## create extinct column - opposite of persisted so the values from the model are probability of extinction, not persistance
## 1 = species is locally extinct at that site
## 0 = species has persisted at that site
nmrsdata_well_expand2$`2012-2016` <- ifelse(nmrsdata_well_expand2$`2012-2016`==0, 1, 0) 
colnames(nmrsdata_well_expand2)[4] <- "extinct" # change column name to extinct
nmrsdata_well_expand2[3] <- NULL # remove TP1 presence/absence
## save file
write.csv(nmrsdata_well_expand2, file="Data/nmrsdata_well_local_extinction.csv", row.names=FALSE)


####################################
##  heavily recorded hectads only ##
####################################

rm(list = ls())

library(ggplot2)
library(ggeffects)
library(DHARMa)
library(lme4)
library(MuMIn)
library(stringr)
library(blmeco)
library(MuMIn)

## Read in data
# NMRS data for cool-adapted moths with elevation and temperature data at a 10km scale
nmrsdata <- readRDS("Data/NMRS_cool_moths_final.rds") ## NMRS data for all hectads and all years with elevation
## hectad recording levels
hec_records <- read.csv("Data/Hectad_recording_levels_1975_1991_2012_2016.csv", header=TRUE)
# recording effort
recording_effort <- read.csv("Data/Recording_effort.csv", header=TRUE)

nmrsdata <- nmrsdata[!nmrsdata$Common_name=="Lunar Thorn",]

## Step 1a. Filter NMRS data to heavily recorded hectads only
heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
heavy_hecs <- heavy_hecs[which(heavy_hecs$n_row > 1), ]
heavy_hecs$n_row <- NULL ## 1084 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_heavy <- nmrsdata[which(nmrsdata$Hectad %in% heavy_hecs$HECTAD), ]
length(unique(nmrsdata_heavy$Hectad)) ## 592 hectads
length(unique(nmrsdata_heavy$Common_name)) ## 72 species

## Step 1b. Define colonised, persisted, and extirpated hectads
nmrsdata_heavy$Recorded <- 1
nmrsdata_heavy_expand <- nmrsdata_heavy %>% tidyr::expand(Common_name, Time_period, Hectad)
recorded <- nmrsdata_heavy[,c("Common_name", "Time_period", "Hectad", "Recorded")]
nmrsdata_heavy_expand <- merge(nmrsdata_heavy_expand, recorded, by=c("Common_name", "Time_period", "Hectad"), all.x=TRUE)
nmrsdata_heavy_expand[is.na(nmrsdata_heavy_expand)] <- 0 ## change NAs to zero == species was NOT recorded at this hectad in this time period
## change to long format
nmrsdata_heavy_expand2 <- nmrsdata_heavy_expand %>%
  spread(Time_period, Recorded) ## each species has 1424 rows = the number of recorded hectads
## first remove rows where TP1 AND TP2 == 0 (this a hectad where a species was never recorded)
nmrsdata_heavy_expand2<-nmrsdata_heavy_expand2[!(nmrsdata_heavy_expand2$`1975-1991`==0 & nmrsdata_heavy_expand2$`2012-2016`==0),]
nmrsdata_heavy_expand2$Hectad_category <- case_when(
  nmrsdata_heavy_expand2$`1975-1991`==0 & nmrsdata_heavy_expand2$`2012-2016`==1 ~ "Colonisation",
  nmrsdata_heavy_expand2$`1975-1991`==1 & nmrsdata_heavy_expand2$`2012-2016`==0 ~ "Extirpation",
  TRUE ~ "Persistence"
)

## Step 1c. Plot the proportion of extinctions at each hectad on a map
## add lat/lon data back in
lat_lon <- nmrsdata_heavy[,c("Hectad", "lat", "lon","elevation10x10km","elevation10x10km_SD")]
lat_lon <- lat_lon %>% distinct(Hectad, .keep_all = TRUE) ## different way of unique
nmrsdata_heavy_expand2 <- merge(nmrsdata_heavy_expand2, lat_lon, by="Hectad", all.x=TRUE)

## find proportion of species at each extirpated hectad - where are extirpations occurring? 
extirpated_hecs1 <- nmrsdata_heavy_expand2 %>% group_by(Hectad, lat, lon) %>% filter(`1975-1991`==1) %>%
  dplyr::summarise(tot_sp = n()) # number of species recorded in each hectad in TP1
extirpated_hecs2 <- nmrsdata_heavy_expand2 %>% group_by(Hectad, lat, lon) %>% filter(sum(`1975-1991`)==0) %>%
  dplyr::summarise(tot_sp = 0)
extirpated_hecs3 <- nmrsdata_heavy_expand2 %>% group_by(Hectad, lat, lon) %>%
  filter(Hectad_category=="Extirpation") %>%
  dplyr::summarise(extir_sp = n()) # number of species which have gone extinct from each hectad
extirpated_hecs <- rbind(extirpated_hecs1, extirpated_hecs2)
extirpated_hecs <- merge(extirpated_hecs, extirpated_hecs3, by=c("Hectad", "lat", "lon"), all.x=TRUE)
extirpated_hecs[is.na(extirpated_hecs)] <- 0 ## change NAs to zero - where there are no extirpations at a hectad
extirpated_hecs$extir_prop <- extirpated_hecs$extir_sp/extirpated_hecs$tot_sp ## 592
extirpated_hecs$extir_prop[is.nan(extirpated_hecs$extir_prop)]<-0 # convert NaNs to zero (hectads where there are only colonising species - the 45 hectads in extirpated_hecs2)

## map
worldmap = map_data('world')
heavy_extir_hecs <- ggplot() +
  geom_polygon(data = worldmap,
               aes(x = long, y = lat, group = group),
               fill = 'gray90', color = 'black') +
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) +
  theme_void() +
  geom_point(data =extirpated_hecs,
             aes(x = as.numeric(lon),
                 y = as.numeric(lat), colour=extir_prop), size=1) +
  scale_color_viridis_c(name="Proportion of species") + 
  theme(title = element_text(size = 12)) +
  scale_fill_continuous(guide = guide_legend()) +
  guides(colour = guide_colourbar(title.vjust = 0.9)) +
  theme(legend.position="bottom", legend.key.height = unit(0.5, "cm"), legend.key.width = unit(2, "cm"))
heavy_extir_hecs
ggsave(heavy_extir_hecs, file="Outputs/Maps/Heavy_extirpated_hecs.png")

## put maps together
legend <- cowplot::get_legend(well_extir_hecs + theme(legend.position = "bottom"))
library(gridExtra)
library(grid)
lay <- rbind(c(1,2), c(3,3))
th <- sum(legend$heights)
extir_maps <- grid.arrange(well_extir_hecs + theme(legend.position="none"), 
                           heavy_extir_hecs + theme(legend.position="none"), legend, 
                           ncol=2, layout_matrix=lay, heights = unit.c(unit(1, "null"), th))
extir_maps
ggsave("Outputs/Maps/Extir_maps.png", extir_maps, height=6, width=10)

## Step 1d. Add in recording effort and climate data and calculate difference between T1 and T2 temperature and precipitation
# add recording effort in
nmrsdata_heavy_expand2 <- merge(nmrsdata_heavy_expand2, recording_effort, by="Hectad", all.x=TRUE)
## add in climate data to each data frame
nmrs_climate <- readRDS("Data/NMRS_hectad_elevation_climate.rds")
climate_hecs <- unique(nmrs_climate[,c(1,2,16,17)]) # 2719 rows
# change from long to wide
climate_wide <- reshape(climate_hecs, idvar = "Hectad", timevar = "Time_period", direction = "wide")
heavy_recs <- unique(nmrsdata_heavy$Hectad)
climate_wide2 <- climate_wide[which(climate_wide$Hectad %in% heavy_recs), ]
colnames(climate_wide2)[2:5] <- c("temperature.TP1", "total_precip.TP1", "temperature.TP2", "total_precip.TP2")
climate_wide2$temp_diff <- climate_wide2$temperature.TP2 - climate_wide2$temperature.TP1
climate_wide2$precip_diff <- climate_wide2$total_precip.TP2 - climate_wide2$total_precip.TP1
## merge in with nmrs data
nmrsdata_heavy_expand2 <- merge(nmrsdata_heavy_expand2, climate_wide2, by="Hectad", all.x=TRUE)

## Step 1e. Remove colonisations - only interested in extinctions and persistence 
nmrsdata_heavy_expand2<-nmrsdata_heavy_expand2[!(nmrsdata_heavy_expand2$Hectad_category=="Colonisation"),] ## remove colonisations
## create extinct column - opposite of TP2 presence/absence so the values from the model are probability of extinction, not persistance
## 1 = species is locally extinct at that site
## 0 = species has persisted at that site
nmrsdata_heavy_expand2$`2012-2016` <- ifelse(nmrsdata_heavy_expand2$`2012-2016`==0, 1, 0) 
colnames(nmrsdata_heavy_expand2)[4] <- "extinct" # change column name to extinct
nmrsdata_heavy_expand2[3] <- NULL # remove TP1 presence/absence
## save file
write.csv(nmrsdata_heavy_expand2, file="Data/nmrsdata_heavy_local_extinction.csv", row.names=FALSE)


#########################################################################################################
#########################################################################################################
#########################################################################################################

########################## Step 2: Run binomial GLMMs ########################## 

#####################################
## well & heavily recorded hectads ##
#####################################

nmrs_glmm <- read.csv("Data/nmrsdata_well_local_extinction.csv", header=TRUE)

# check correlations between climate variables
round(cor(nmrs_glmm[,c("temperature.TP1","temp_diff","total_precip.TP1","precip_diff", "n_recs")]),3)
# all <0.7
## leave all variables in
round(cor(nmrs_glmm[,c("temperature.TP2","temp_diff","total_precip.TP2","precip_diff", "n_recs")]),3)
# 0.62 temp_diff and temp_tp2
# 0.66 precip_diff and precip_tp2
## leave all in

str(nmrs_glmm)

## TP1 model 
model_tp1 <- glmer(extinct ~ scale(temperature.TP1)*scale(total_precip.TP1) + scale(temp_diff)*scale(precip_diff) +
                 scale(temperature.TP1)*scale(precip_diff) + scale(total_precip.TP1)*scale(temp_diff) +
                 scale(total_precip.TP1)*scale(precip_diff) + scale(temperature.TP1)*scale(temp_diff) +
                 scale(n_recs) + (1|Common_name), family="binomial", data=nmrs_glmm, 
                 na.action = "na.fail")
AIC(model_tp1)
summary(model_tp1) # n_records highly significant (negative)
r.squaredGLMM(model_tp1) 
car::vif(model_tp1) ## all < 3

## check model assumptions - only use this to check for outliers, binomial GLMM don't require normal residuals or homogeneity of variance
testDispersion(model_tp1) ## underdispersed (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = model_tp1, plot = F)
plot(simulationOutput) ## no outliers

overdisp_fun <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  } 
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model, type = "pearson") # computes pearson residuals
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
} 
overdisp_fun(model_tp1) # 0.98 no issues
dispersion_glmer(model_tp1) # 1.1 - suggests not an issue

## test for spatial autocorrelation in residuals of model
# first get unique coordinates
nmrs_glmm$coords <- paste(nmrs_glmm$lon,", ",nmrs_glmm$lat)
coords <- c(unique(nmrs_glmm$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))
sims<-simulateResiduals(model_tp1)
simsrecalc<-recalculateResiduals(sims,group = nmrs_glmm$coords) # recalculate residuals based on unique coordinates
testSpatialAutocorrelation(simsrecalc, x = x_unique, y = y_unique) ## non significant - no evidence of spatial autocorrelation

## dredge model
# Dredge & get best models based on AIC
## dredge & r2 doesn't work with an optimizer 
## try this instead (from: https://stackoverflow.com/questions/53856379/dredge-doesnt-work-when-specifying-glmer-optimizer)

# Fit a null model with RE (use a non-exported function .nullFitRE):
nullmodel <- MuMIn:::.nullFitRE(model_tp1)
# the above step is not necessary, but avoids repeated re-fitting of the null model.

d1_tp1 <- dredge(model_tp1)

## take models with change in AICc <6 and average them
topmods_tp1 <- subset(d1_tp1, delta <6)
topmods_tp1_df <- data.frame(topmods_tp1)
# save 
write.csv(topmods_tp1_df, file="Outputs/Results/Candidate_models_TP1_well_receffort.csv", row.names=FALSE)
avgmod_tp1 <- model.avg(d1_tp1, fit=TRUE) 
summary(avgmod_tp1)


## save summary
avg_tp1_summary <- summary(avgmod_tp1, full=F) #pulling out model averages
tp1_summary<-as.data.frame(avg_tp1_summary$coefmat.full) #selecting full model coefficient averages
## get confidence intervals
confint(avgmod_tp1, full=F)
CI <- as.data.frame(confint(avgmod_tp1, full=F)) # get confidence intervals for full model
tp1_summary$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
tp1_summary$CI.max <-CI$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(tp1_summary, keep.rownames = "coefficient") #put rownames into column
names(tp1_summary) <- gsub(" ", "", names(tp1_summary)) # remove spaces from column headers
## get importance values
importance_tp1 <- MuMIn::sw(avgmod_tp1) ## total precip, temperature and temperature difference in all 34 models
importance_tp1 <- as.data.frame(importance_tp1)
importance_tp1$parameters <- row.names(importance_tp1)
row.names(importance_tp1) <- 1:nrow(importance_tp1)
tp1_summary <- merge(tp1_summary, importance_tp1, by.x="coefficient", by.y="parameters", all=TRUE)
write.csv(tp1_summary, file="Outputs/Results/Average_model_summary_TP1_well_receffort.csv", row.names=FALSE)

nmrs_glmm <- nmrs_glmm[,c("Hectad", "Common_name", "temperature.TP1", "total_precip.TP1", "temp_diff", "precip_diff", "extinct")]

## Heatmap plots ##
# temp TP1 * precip TP1
pred.dat_tp1 = expand.grid(total_precip.TP1 = seq(min(nmrs_glmm$total_precip.TP1), max(nmrs_glmm$total_precip.TP1), length=100), 
                       temperature.TP1 =  seq(min(nmrs_glmm$temperature.TP1), max(nmrs_glmm$temperature.TP1), length=100),
                       precip_diff = median(nmrs_glmm$precip_diff),temp_diff = median(nmrs_glmm$temp_diff),
                       Common_name=unique(nmrs_glmm$Common_name))
# Add predictions
pred.dat_tp1$extinction = predict(avgmod_tp1, newdata=pred.dat_tp1, type="response", re.form=NULL)
write.csv(pred.dat_tp1, file="Outputs/Results/Predicted_extinction_tp1_well_sp.csv", row.names=FALSE)
library(insight)
library(ggeffects)

install.packages("insight")
devtools::install_github("strengejacke/ggeffects")

pred2 <- ggpredict(avgmod_tp1, terms=c("temperature.TP1", "total_precip.TP1"), type="re")

for(i in unique(pred.dat_tp1$Common_name)) {
  print(i)
  # Printing ggplot within for-loop
  
  #   theme_classic()
  nmrs_sp <- nmrs_glmm[nmrs_glmm$Common_name==i,]
  nmrs_sp2 <- unique(nmrs_sp[,c("Hectad","total_precip.TP1", "temperature.TP1")])
  
  extinct_plot <- ggplot(pred.dat_tp1[pred.dat_tp1$Common_name==i,], aes(total_precip.TP1, temperature.TP1, fill=extinction)) +
    geom_tile() +
    scale_fill_gradient2(low="dodgerblue2", mid="white", high="red",
                         midpoint=0.5) +
    labs(fill="Probability of \nlocal extinction") +
    geom_point(data=nmrs_sp2, aes(total_precip.TP1, temperature.TP1), colour="grey", alpha=0.5, inherit.aes = FALSE) +
    xlab(expression("Total precipitation T"[1])) + 
    ylab(expression("Mean temperature T"[1])) + 
    ggtitle(i) +
    coord_cartesian(xlim = c(quantile(nmrs_glmm$total_precip.TP1, 0.05), quantile(nmrs_glmm$total_precip.TP1, 0.95)),
                    ylim = c(quantile(nmrs_glmm$temperature.TP1, 0.05), quantile(nmrs_glmm$temperature.TP1, 0.95))) +
    theme_classic()
  
  ggsave(extinct_plot, file=paste0("Outputs/Graphs/Extinction_TP1/Extinction_prob_TP1_", i,".png"), width = 20, height = 15, units = "cm")
  Sys.sleep(2)
}


## Population level predictions
pred.dat_tp1_pop = expand.grid(total_precip.TP1 = seq(min(nmrs_glmm$total_precip.TP1), max(nmrs_glmm$total_precip.TP1), length=100), 
                               temperature.TP1 =  seq(min(nmrs_glmm$temperature.TP1), max(nmrs_glmm$temperature.TP1), length=100),
                               precip_diff = median(nmrs_glmm$precip_diff),temp_diff = median(nmrs_glmm$temp_diff), 
                               n_recs=median(nmrs_glmm$n_recs), Common_name=nmrs_glmm$Common_name[1])

# Add predictions
pred.dat_tp1_pop$extinction = predict(avgmod_tp1, newdata=pred.dat_tp1_pop, type="response", re.form=NA) ## mean population predictions
write.csv(pred.dat_tp1_pop, file="Outputs/Results/Predicted_extinction_tp1_well_pop.csv", row.names=FALSE)

pred.dat_tp1_pop <- read.csv("Outputs/Results/Predicted_extinction_tp1_well_pop.csv", header=TRUE)

nmrs_temp_precip <- unique(nmrs_glmm[,c("Hectad","total_precip.TP1", "temperature.TP1")])
precip_temp_TP1 <- ggplot(pred.dat_tp1_pop, aes(total_precip.TP1, temperature.TP1, fill=extinction)) +
  geom_tile() +
  scale_fill_gradient2(low="dodgerblue2", mid="white", high="red",
                       midpoint=0.5, limits=c(0,1),
                       labels = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  labs(fill="Probability of \nlocal extinction") +
  geom_point(data=nmrs_temp_precip, aes(total_precip.TP1, temperature.TP1), colour="grey", alpha=0.5, inherit.aes = FALSE) +
  xlab(expression("Total precipitation T"[1])) + 
  ylab(expression("Mean temperature T"[1])) + 
  coord_cartesian(xlim = c(quantile(nmrs_glmm$total_precip.TP1, 0.05), quantile(nmrs_glmm$total_precip.TP1, 0.95)),
                  ylim = c(quantile(nmrs_glmm$temperature.TP1, 0.05), quantile(nmrs_glmm$temperature.TP1, 0.95))) +
  theme_classic() +
  theme(legend.key.height = unit(1.2, "cm"), legend.key.width = unit(0.8, "cm")) 
## 1 = persist, 0 = extirpate
## red = high probability of extirpatation
## blue = high probability of persistance
precip_temp_TP1
ggsave(precip_temp_TP1, file="Outputs/Graphs/Temp_precip_TP1_heatmap_well_receffort.png", height=5, width=7)



#############################################################################
## TP2 model 
model_tp2 <- glmer(extinct ~ scale(temperature.TP2) + scale(temp_diff) + scale(total_precip.TP2) + scale(precip_diff) + 
                 scale(temperature.TP2)*scale(total_precip.TP2) + scale(temp_diff)*scale(precip_diff) +
                 scale(temperature.TP2)*scale(precip_diff) + scale(total_precip.TP2)*scale(temp_diff) +
                 scale(total_precip.TP2)*scale(precip_diff) + scale(temperature.TP2)*scale(temp_diff) +
                 scale(n_recs) + (1|Common_name), family="binomial", glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                 data=nmrs_glmm, na.action = "na.fail")
summary(model_tp2)

library(MuMIn)
r.squaredGLMM(model_tp2) 
car::vif(model_tp2) ## all < 4

## check model assumptions - only use this to check for outliers, binomial GLMM don't require normal residuals or homogeneity of variance
testDispersion(model_tp2) ## slight overdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = model_tp2, plot = F)
plot(simulationOutput) ## no outliers

library(blmeco)
dispersion_glmer(model_tp2) # 1.1 - suggests not an issue
overdisp_fun <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  } 
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model, type = "pearson") # computes pearson residuals
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
} 
overdisp_fun(model_tp2) # 0.98

## test for spatial autocorrelation in residuals of model
# first get unique coordinates
library(stringr)
nmrs_glmm$coords <- paste(nmrs_glmm$lon,", ",nmrs_glmm$lat)
coords <- c(unique(nmrs_glmm$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))
sims<-simulateResiduals(model_tp2)
simsrecalc<-recalculateResiduals(sims,group = nmrs_glmm$coords) # recalculate residuals based on unique coordinates
testSpatialAutocorrelation(simsrecalc, x = x_unique, y = y_unique) ## non significant - no evidence of spatial autocorrelation

## dredge model
# Dredge & get best models based on AIC
library(MuMIn)
## dredge & r2 doesn't work with an optimizer 
## try this instead (from: https://stackoverflow.com/questions/53856379/dredge-doesnt-work-when-specifying-glmer-optimizer)

# Fit a null model with RE (use a non-exported function .nullFitRE):
nullmodel <- MuMIn:::.nullFitRE(model_tp2)
# the above step is not necessary, but avoids repeated re-fitting of the null model.

d1_tp2 <- dredge(model_tp2, beta="partial.sd", extra =list(R2 = function(x) {
  r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

## take models with change in AICc <6 and average them
topmods_tp2 <- subset(d1_tp2, delta <6)
topmods_tp2_df <- data.frame(topmods_tp2)
# save 
write.csv(topmods_tp2_df, file="Outputs/Results/Candidate_models_TP2_well_receffort.csv", row.names=FALSE)
avgmod_tp2 <- model.avg(topmods_tp2, fit=TRUE) 
summary(avgmod_tp2)

## save summary
avg_tp2_summary <- summary(avgmod_tp2, full=F) #pulling out model averages
tp2_summary<-as.data.frame(avg_tp2_summary$coefmat.full) #selecting full model coefficient averages
## get confidence intervals
confint(avgmod_tp2, full=F)
CI <- as.data.frame(confint(avgmod_tp2, full=F)) # get confidence intervals for full model
tp2_summary$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
tp2_summary$CI.max <-CI$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(tp2_summary, keep.rownames = "coefficient") #put rownames into column
names(tp2_summary) <- gsub(" ", "", names(tp2_summary)) # remove spaces from column headers
## get importance values
tp2_summary$parameters <- row.names(tp2_summary)
row.names(tp2_summary) <- 1:nrow(tp2_summary)
importance_tp2 <- MuMIn::sw(avgmod_tp2) ## total precip, temperature and temperature difference in all 34 models
importance_tp2 <- as.data.frame(importance_tp2)
importance_tp2$parameters <- row.names(importance_tp2)
row.names(importance_tp2) <- 1:nrow(importance_tp2)
tp2_summary <- merge(tp2_summary, importance_tp2, by.x="coefficient", by.y="parameters", all=TRUE)
write.csv(tp2_summary, file="Outputs/Results/Average_model_summary_TP2_well_receffort.csv", row.names=FALSE)

## Heatmap plots ##
# Interaction between annual temperature & total precipitation
# First for each species and plot heatmaps to see if mean population trend is similar across all species
# Remember that species which only occur at a few sites will be very close to pop mean as cannot tell what is happening from few data points
# Model takes mean first, then adjusts for each site - so few sites = closer to pop mean
pred.dat_tp2 = expand.grid(temp_diff = median(nmrs_glmm$temp_diff), precip_diff = median(nmrs_glmm$precip_diff),
                       temperature.TP2 = seq(min(nmrs_glmm$temperature.TP2), max(nmrs_glmm$temperature.TP2), length=100),
                       total_precip.TP2 = seq(min(nmrs_glmm$total_precip.TP2), max(nmrs_glmm$total_precip.TP2), length=100),
                       Common_name=unique(nmrs_glmm$Common_name))

# Add predictions
pred.dat_tp2$extinction = predict(avgmod_tp2, newdata=pred.dat_tp2, type="response", re.form=NULL)
write.csv(pred.dat_tp2, file="Outputs/Results/Predicted_extinction_tp2_well_sp.csv", row.names=FALSE)

for(i in unique(pred.dat_tp2$Common_name)) {
  print(i)
  # Printing ggplot within for-loop

  #   theme_classic()
  nmrs_sp <- nmrs_glmm[nmrs_glmm$Common_name==i,]
  nmrs_sp2 <- unique(nmrs_sp[,c("Hectad","total_precip.TP2", "temperature.TP2")])
  
  extinct_plot <- ggplot(pred.dat_tp2[pred.dat_tp2$Common_name==i,], aes(total_precip.TP2, temperature.TP2, fill=extinction)) +
    geom_tile() +
    scale_fill_gradient2(low="dodgerblue2", mid="white", high="red",
                         midpoint=0.5) +
    labs(fill="Probability of \nlocal extinction") +
    geom_point(data=nmrs_sp2, aes(total_precip.TP2, temperature.TP2), colour="grey", alpha=0.5, inherit.aes = FALSE) +
    xlab(expression("Total precipitation T"[2])) + 
    ylab(expression("Mean temperature T"[2])) + 
    ggtitle(i) +
    coord_cartesian(xlim = c(quantile(nmrs_sp2$total_precip.TP2, 0.05), quantile(nmrs_sp2$total_precip.TP2, 0.95)),
                    ylim = c(quantile(nmrs_sp2$temperature.TP2, 0.05), quantile(nmrs_sp2$temperature.TP2, 0.95))) +
    theme_classic()
  
  ggsave(extinct_plot, file=paste0("Outputs/Graphs/Extinction_TP2_2/Extinction_prob_TP2_", i,".png"), width = 20, height = 15, units = "cm")
  Sys.sleep(2)
}

## Population level predictions
pred.dat_tp2_pop = expand.grid(temp_diff = median(nmrs_glmm$temp_diff), precip_diff = median(nmrs_glmm$precip_diff), 
                           temperature.TP2 = seq(min(nmrs_glmm$temperature.TP2), max(nmrs_glmm$temperature.TP2), length=100),
                           total_precip.TP2 = seq(min(nmrs_glmm$total_precip.TP2), max(nmrs_glmm$total_precip.TP2), length=100),
                           Common_name=nmrs_glmm$Common_name[1], n_recs=median(nmrs_glmm$n_recs))

# Add predictions
pred.dat_tp2_pop$extinction = predict(avgmod_tp2, newdata=pred.dat_tp2_pop, type="response", re.form=NA)
write.csv(pred.dat_tp2_pop, file="Outputs/Results/Predicted_extinction_tp2_well_pop.csv", row.names=FALSE)

pred.dat_tp2_pop <- read.csv("Outputs/Results/Predicted_extinction_tp2_well_pop.csv", header=TRUE)

nmrs_temp_precip <- unique(nmrs_glmm[,c("Hectad","total_precip.TP2", "temperature.TP2")])
precip_temp_TP2 <- ggplot(pred.dat_tp2_pop, aes(total_precip.TP2, temperature.TP2, fill=extinction)) +
  geom_tile() +
  scale_fill_gradient2(low="dodgerblue2", mid="white", high="red",
                       midpoint=0.5, limits=c(0,1),
                       labels = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  labs(fill="Probability of \nlocal extinction") +
  geom_point(data=nmrs_temp_precip, aes(total_precip.TP2, temperature.TP2), colour="grey", alpha=0.5, inherit.aes = FALSE) +
  xlab(expression("Total precipitation T"[2])) + 
  ylab(expression("Mean temperature T"[2])) + 
  coord_cartesian(xlim = c(quantile(nmrs_glmm$total_precip.TP2, 0.05), quantile(nmrs_glmm$total_precip.TP2, 0.95)),
                  ylim = c(quantile(nmrs_glmm$temperature.TP2, 0.05), quantile(nmrs_glmm$temperature.TP2, 0.95))) +
  theme_classic() +
  theme(legend.key.height = unit(1.2, "cm"), legend.key.width = unit(0.8, "cm")) 
## 1 = persist, 0 = extirpate
## red = high probability of extirpatation
## blue = high probability of persistance
precip_temp_TP2
ggsave(precip_temp_TP2, file="Outputs/Graphs/Temp_precip_TP2_heatmap_well.png", height=5, width=7)

## calculate some extinction probabilities for paper

temp95 <- quantile(nmrs_glmm$temperature.TP2, 0.95)
mean_temp <- mean(nmrs_glmm$temperature.TP2)
precip95 <- quantile(nmrs_glmm$total_precip.TP2, 0.95)
precip05 <- quantile(nmrs_glmm$total_precip.TP2, 0.05)

pred.dat = expand.grid(temp_diff = median(nmrs_glmm$temp_diff), precip_diff = median(nmrs_glmm$precip_diff), 
                               temperature.TP2 = c(mean_temp,mean_temp, temp95,temp95),
                               total_precip.TP2 = c(precip05, precip95,precip05, precip95),
                               Common_name=nmrs_glmm$Common_name[1], n_recs=median(nmrs_glmm$n_recs))
pred.dat <- unique(pred.dat)
pred.dat$extinction <-predict(avgmod_tp2, newdata = pred.dat, type="response", re.form=NA)

tp1_pop <- pred.dat_tp1_pop[(pred.dat_tp1_pop$temperature.TP1 >= temp05 & pred.dat_tp1_pop$temperature.TP1 <= temp95), ]
tp1_pop <- tp1_pop[(tp1_pop$total_precip.TP1 >= precip05 & tp1_pop$total_precip.TP1 <= precip95), ]

## Forest plots for main effects (TP2 only)
tp2_summary <- read.csv("Outputs/Results/Average_model_summary_TP2_well_receffort.csv", header=TRUE)

tp2_summary <- tp2_summary[-1,]
tp2_summary <- tp2_summary[,c(1:2,7:8)]

lookup <- c("scale(precip_diff)"="\U0394 Precipitation", "scale(temp_diff)"="\U0394 Temperature",
            "scale(temperature.TP1)"="Temperature T1","scale(total_precip.TP1)"="Precipitation T1",
            "scale(temperature.TP2)"="Temperature T2", "scale(total_precip.TP2)"="Precipitation T2",
            "scale(precip_diff):scale(temp_diff)"="\U0394 Precipitation x \U0394 Temperature",
            "scale(precip_diff):scale(temperature.TP1)"="\U0394 Precipitation x Temperature T1",
            "scale(precip_diff):scale(total_precip.TP1)" = "\U0394 Precipitation x Precipitation T1",
            "scale(temp_diff):scale(temperature.TP1)" = "\U0394 Temperature x Temperature T1",
            "scale(temp_diff):scale(total_precip.TP1)"="\U0394 Temperature x Precipitation T1",
            "scale(temperature.TP1):scale(total_precip.TP1)"="Temperature T1 x Precipitation T1",
            "scale(precip_diff):scale(temperature.TP2)" = "\U0394 Precipitation x Temperature T2",
            "scale(precip_diff):scale(total_precip.TP2)"="\U0394 Precipitation x Precipitation T2",
            "scale(temp_diff):scale(temperature.TP2)" = "\U0394 Temperature x Temperature T2",
            "scale(temp_diff):scale(total_precip.TP2)"="\U0394 Temperature x Precipitation T2",
            "scale(temperature.TP2):scale(total_precip.TP2)"="Temperature T2 x Precipitation T2",
            "scale(n_recs)"="Number of records")
            
tp2_summary <- tp2_summary %>% mutate(coefficient=as.list(lookup[as.character(coefficient)]))

row.names(tp2_summary) <- 1:nrow(tp2_summary)
tp2_summary <- tp2_summary[c(9,11,2,6,10,3,4,5,7,8,1),]
tp2_summary$coefficient <- as.character(tp2_summary$coefficient)
tp2_summary$coefficient <- factor(tp2_summary$coefficient, levels=c("Number of records", "\U0394 Temperature x Precipitation T2", "\U0394 Temperature x Temperature T2",
                                                            "\U0394 Precipitation x Precipitation T2", "\U0394 Precipitation x Temperature T2",
                                                            "\U0394 Precipitation x \U0394 Temperature", "Temperature T2 x Precipitation T2",
                                                           "\U0394 Temperature", "\U0394 Precipitation", "Precipitation T2", "Temperature T2"))


tp2_summary <- tp2_summary %>% mutate(colour = case_when(
  CI.min < 0 & CI.max < 0 ~ "yes",
  CI.min > 0 & CI.max > 0 ~ "yes",
  TRUE ~ "no"
))

forest_plot_well_tp2 <- ggplot(data=tp2_summary, aes(x=coefficient, y=Estimate, colour=colour))+ #excluding intercept because estimates so much larger
  geom_point(size=4)+ #points for coefficient estimates
  theme_classic(base_size = 12)+ #clean graph
  geom_errorbar(aes(ymin=CI.min, ymax=CI.max, colour=colour), # CIs
                width=.4,lwd=0.6) +
  geom_hline(yintercept=0, linetype='dotted', col = 'black')+
  labs(y="Scaled model-averaged parameter estimates", x=" ") +
  theme(axis.text=element_text(color="black"), axis.text.y = element_text(size = 12)) +
  scale_colour_manual(values = c("yes" = "black", "no" = "lightgrey")) +
  theme(legend.position="none") +
  coord_flip()
forest_plot_well_tp2
ggsave(forest_plot_well_tp2, file="Outputs/Graphs/Forest_plot_well_tp2.png", height=6, width=7)

## put the heatmap plot and forest plot for TP2 together as Figure 3
legend <- cowplot::get_legend(precip_temp_TP2 + theme(legend.position = "bottom"))
library(gridExtra)
library(grid)
library(ggpubr)
lay <- rbind(c(1,2), c(3,3))
th <- sum(legend$heights)
figure3 <- ggpubr::ggarrange(forest_plot_well_tp2, precip_temp_TP2, labels=c("(a)", "(b)"),
                             nrow=1, ncol=2)
figure3
ggsave("Outputs/Graphs/Figure3_well_tp2_receffort.png", figure3, height=6, width=14)


#######################################################
### Histograms of mean temp and total precipitation at each hectad for supplementary material

nmrs_climate <- nmrs_glmm[,c("Hectad", "temperature.TP1", "temperature.TP2", "total_precip.TP1",
                             "total_precip.TP2", "temp_diff", "precip_diff")]
nmrs_climate <- unique(nmrs_climate) # 869
## take out the lowest and highest 5% of values for each as this is what is plotted on the heatmaps

temp_tp1 <- ggplot(nmrs_climate, aes(x=temperature.TP1)) + geom_histogram() +
  coord_cartesian(xlim = c(quantile(nmrs_climate$temperature.TP1, 0.05), quantile(nmrs_climate$temperature.TP1, 0.95))) +
  xlab(expression("Mean temperature T"[1])) +
  ylab("Number of hectads") +
  #scale_x_continuous(limits = c(5,12)) +
  theme_classic()+
  theme(text = element_text(size = 30))
temp_tp1

precip_tp1 <- ggplot(nmrs_climate, aes(x=total_precip.TP1)) + geom_histogram() +
  coord_cartesian(xlim = c(quantile(nmrs_climate$total_precip.TP1, 0.05), quantile(nmrs_climate$total_precip.TP1, 0.95))) +
  xlab(expression("Total precipitation T"[1])) +
  ylab("Number of hectads") +
  #scale_x_continuous(limits = c(5,12)) +
  theme_classic()+
  theme(text = element_text(size = 30))
precip_tp1

temp_tp2 <- ggplot(nmrs_climate, aes(x=temperature.TP2)) + geom_histogram() +
  coord_cartesian(xlim = c(quantile(nmrs_climate$temperature.TP2, 0.05), quantile(nmrs_climate$temperature.TP2, 0.95))) +
  xlab(expression("Mean temperature T"[2])) +
  ylab("Number of hectads") +
  #scale_x_continuous(limits = c(5,12)) +
  theme_classic()+
  theme(text = element_text(size = 30))
temp_tp2

precip_tp2 <- ggplot(nmrs_climate, aes(x=total_precip.TP2)) + geom_histogram() +
  coord_cartesian(xlim = c(quantile(nmrs_climate$total_precip.TP2, 0.05), quantile(nmrs_climate$total_precip.TP2, 0.95))) +
  xlab(expression("Total precipitation T"[2])) +
  ylab("Number of hectads") +
  #scale_x_continuous(limits = c(5,12)) +
  theme_classic()+
  theme(text = element_text(size = 30))
precip_tp2

# add in temperature change and precipitation change to these
temp_diff <- ggplot(nmrs_climate, aes(x=temp_diff)) + geom_histogram() +
  coord_cartesian(xlim = c(quantile(nmrs_climate$temp_diff, 0.05), quantile(nmrs_climate$temp_diff, 0.95))) +
  xlab("Temperature change") +
  ylab("Number of hectads") +
  #scale_x_continuous(limits = c(5,12)) +
  theme_classic()+
  theme(text = element_text(size = 30))
temp_diff

precip_diff <- ggplot(nmrs_climate, aes(x=precip_diff)) + geom_histogram() +
  coord_cartesian(xlim = c(quantile(nmrs_climate$precip_diff, 0.05), quantile(nmrs_climate$precip_diff, 0.95))) +
  xlab("Precipitation change") +
  ylab("Number of hectads") +
  #scale_x_continuous(limits = c(5,12)) +
  theme_classic()+
  theme(text = element_text(size = 30))
precip_diff

climate_hist <- ggarrange(temp_tp1, precip_tp1, temp_tp2, precip_tp2, temp_diff, precip_diff, 
                          labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), nrow=3, ncol=2,
                          font.label = list(size = 30), vjust=-0.6)+
  theme(plot.margin = margin(1.4,0.1,0.1,0.1, "cm")) 


climate_hist 
ggsave(climate_hist, file="Outputs/Graphs/Climate_histograms.png", height=20, width=15)



##################################################################################################################
##################################################################################################################
##################################################################################################################

###################################
## heavily recorded hectads only ##
###################################

nmrs_glmm <- read.csv("Data/nmrsdata_heavy_local_extinction.csv", header=TRUE)

# check correlations between climate variables
round(cor(nmrs_glmm[,c("temperature.TP1","temp_diff","total_precip.TP1","precip_diff", "n_recs")]),3)
# all <0.7
## leave all variables in
round(cor(nmrs_glmm[,c("temperature.TP2","temp_diff","total_precip.TP2","precip_diff", "n_recs")]),3)
# 0.62 temp_diff and temp_tp2
# 0.66 precip_diff and precip_tp2
## leave all in
str(nmrs_glmm)

## TP1 model 
model_tp1 <- glmer(extinct ~ scale(temperature.TP1) + scale(temp_diff) + scale(total_precip.TP1) + scale(precip_diff) + 
                     scale(temperature.TP1)*scale(total_precip.TP1) + scale(temp_diff)*scale(precip_diff) +
                     scale(temperature.TP1)*scale(precip_diff) + scale(total_precip.TP1)*scale(temp_diff) +
                     scale(total_precip.TP1)*scale(precip_diff) + scale(temperature.TP1)*scale(temp_diff) +
                     scale(n_recs) + (1|Common_name), family="binomial", data=nmrs_glmm, 
                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                   na.action = "na.fail")
summary(model_tp1)
## temperature*temp_diff is significant
r.squaredGLMM(model_tp1) ## 8.3% fixed effects, 19.9% fixed and random effects
car::vif(model_tp1) ## all < 3

## check model assumptions - only use this to check for outliers, binomial GLMM don't require normal residuals or homogeneity of variance
testDispersion(model_tp1) ## slight overdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = model_tp1, plot = F)
plot(simulationOutput) ## no outliers

overdisp_fun <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  } 
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model, type = "pearson") # computes pearson residuals
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
} 
overdisp_fun(model_tp1) # 0.98 no issues
dispersion_glmer(model_tp1) # 1.1 - suggests not an issue

## test for spatial autocorrelation in residuals of model
# first get unique coordinates
nmrs_glmm$coords <- paste(nmrs_glmm$lon,", ",nmrs_glmm$lat)
coords <- c(unique(nmrs_glmm$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))
sims<-simulateResiduals(model_tp1)
simsrecalc<-recalculateResiduals(sims,group = nmrs_glmm$coords) # recalculate residuals based on unique coordinates
testSpatialAutocorrelation(simsrecalc, x = x_unique, y = y_unique) ## non significant - no evidence of spatial autocorrelation

## dredge model
# Dredge & get best models based on AIC
## dredge & r2 doesn't work with an optimizer 
## try this instead (from: https://stackoverflow.com/questions/53856379/dredge-doesnt-work-when-specifying-glmer-optimizer)

# Fit a null model with RE (use a non-exported function .nullFitRE):
nullmodel <- MuMIn:::.nullFitRE(model_tp1)
# the above step is not necessary, but avoids repeated re-fitting of the null model.

d1_tp1 <- dredge(model_tp1, beta="partial.sd", extra =list(R2 = function(x) {
  r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

## take models with change in AICc <6 and average them
topmods_tp1 <- subset(d1_tp1, delta <6)
topmods_tp1_df <- data.frame(topmods_tp1)
# save 
write.csv(topmods_tp1_df, file="Outputs/Results/Candidate_models_TP1_heavy_receffort.csv", row.names=FALSE)
avgmod_tp1 <- model.avg(topmods_tp1) 
summary(avgmod_tp1)

## save summary
avg_tp1_summary <- summary(avgmod_tp1, full=F) #pulling out model averages
tp1_summary<-as.data.frame(avg_tp1_summary$coefmat.full) #selecting full model coefficient averages
## get confidence intervals
confint(avgmod_tp1, full=F)
CI <- as.data.frame(confint(avgmod_tp1, full=F)) # get confidence intervals for full model
tp1_summary$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
tp1_summary$CI.max <-CI$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(tp1_summary, keep.rownames = "coefficient") #put rownames into column
names(tp1_summary) <- gsub(" ", "", names(tp1_summary)) # remove spaces from column headers
## get importance values
importance_tp1 <- MuMIn::sw(avgmod_tp1) ## total precip, temperature and temperature difference in all 34 models
importance_tp1 <- as.data.frame(importance_tp1)
importance_tp1$parameters <- row.names(importance_tp1)
row.names(importance_tp1) <- 1:nrow(importance_tp1)
tp1_summary <- merge(tp1_summary, importance_tp1, by.x="coefficient", by.y="parameters", all=TRUE)
write.csv(tp1_summary, file="Outputs/Results/Average_model_summary_TP1_heavy_receffort.csv", row.names=FALSE)

## Heatmap plots ##
# temp TP1 * precip TP1
pred.dat_tp1 = expand.grid(total_precip.TP1 = seq(min(nmrs_glmm$total_precip.TP1), max(nmrs_glmm$total_precip.TP1), length=100), 
                           temperature.TP1 =  seq(min(nmrs_glmm$temperature.TP1), max(nmrs_glmm$temperature.TP1), length=100),
                           precip_diff = median(nmrs_glmm$precip_diff),temp_diff = median(nmrs_glmm$temp_diff),
                           Common_name=unique(nmrs_glmm$Common_name))


# Add predictions
pred.dat_tp1$extinction = predict(avgmod_tp1, newdata=pred.dat_tp1, type="response")
write.csv(pred.dat_tp1, file="Outputs/Results/Predicted_extinction_tp1_heavy.csv", row.names=FALSE)

nmrs_temp_precip <- unique(nmrs_glmm[,c("Hectad","total_precip.TP1", "temperature.TP1")])
precip_temp_TP1 <- ggplot(pred.dat_tp1, aes(total_precip.TP1, temperature.TP1, fill=extinction)) +
  geom_tile() +
  scale_fill_gradient2(low="dodgerblue2", mid="white", high="red",
                       midpoint=0.5) +
  labs(fill="Probability of \nlocal extinction") +
  geom_point(data=nmrs_temp_precip, aes(total_precip.TP1, temperature.TP1), colour="grey", alpha=0.5, inherit.aes = FALSE) +
  xlab(expression("Total precipitation T"[1])) + 
  ylab(expression("Mean temperature T"[1])) + 
  coord_cartesian(xlim = c(quantile(nmrs_glmm$total_precip.TP1, 0.05), quantile(nmrs_glmm$total_precip.TP1, 0.95)),
                  ylim = c(quantile(nmrs_glmm$temperature.TP1, 0.05), quantile(nmrs_glmm$temperature.TP1, 0.95))) +
  theme_classic() +
  theme(legend.key.height = unit(0.5, "cm"), legend.key.width = unit(3, "cm"))
## 1 = persist, 0 = extirpate
## red = high probability of extirpatation
## blue = high probability of persistance
precip_temp_TP1
ggsave(precip_temp_TP1, file="Outputs/Graphs/Temp_precip_TP1_heatmap_heavy.png", height=5, width=7)


#############################################################################
## TP2 model 
model_tp2 <- glmer(extinct ~ scale(temperature.TP2) + scale(temp_diff) + scale(total_precip.TP2) + scale(precip_diff) + 
                     scale(temperature.TP2)*scale(total_precip.TP2) + scale(temp_diff)*scale(precip_diff) +
                     scale(temperature.TP2)*scale(precip_diff) + scale(total_precip.TP2)*scale(temp_diff) +
                     scale(total_precip.TP2)*scale(precip_diff) + scale(temperature.TP2)*scale(temp_diff) +
                     scale(n_recs) + (1|Common_name), family="binomial", glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                   data=nmrs_glmm, na.action = "na.fail")
summary(model_tp2)
## 3 way interaction: total_precip*temp_diff*temperature is non-significant
## temperature*temp_diff is significant

library(MuMIn)
r.squaredGLMM(model_tp2) ## 8.3% fixed effects, 19.8% fixed and random effects
car::vif(model_tp2) ## all < 4

## check model assumptions - only use this to check for outliers, binomial GLMM don't require normal residuals or homogeneity of variance
testDispersion(model_tp2) ## slight overdispersion (if red line is to the left = underdispersion)
simulationOutput <- simulateResiduals(fittedModel = model_tp2, plot = F)
plot(simulationOutput) ## no outliers

library(blmeco)
dispersion_glmer(model_tp2) # 1.1 - suggests not an issue
overdisp_fun <- function(model) {
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  } 
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model, type = "pearson") # computes pearson residuals
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
} 
overdisp_fun(model_tp2) # 0.98

## test for spatial autocorrelation in residuals of model
# first get unique coordinates
library(stringr)
nmrs_glmm$coords <- paste(nmrs_glmm$lon,", ",nmrs_glmm$lat)
coords <- c(unique(nmrs_glmm$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))
sims<-simulateResiduals(model_tp2)
simsrecalc<-recalculateResiduals(sims,group = nmrs_glmm$coords) # recalculate residuals based on unique coordinates
testSpatialAutocorrelation(simsrecalc, x = x_unique, y = y_unique) ## non significant - no evidence of spatial autocorrelation


## dredge model
# Dredge & get best models based on AIC
library(MuMIn)
## dredge & r2 doesn't work with an optimizer 
## try this instead (from: https://stackoverflow.com/questions/53856379/dredge-doesnt-work-when-specifying-glmer-optimizer)

# Fit a null model with RE (use a non-exported function .nullFitRE):
nullmodel <- MuMIn:::.nullFitRE(model_tp2)
# the above step is not necessary, but avoids repeated re-fitting of the null model.

d1_tp2 <- dredge(model_tp2, beta="partial.sd", extra =list(R2 = function(x) {
  r.squaredGLMM(x, null = nullmodel)["delta", ]
}))

## take models with change in AICc <6 and average them
topmods_tp2 <- subset(d1_tp2, delta <6)
topmods_tp2_df <- data.frame(topmods_tp2)
# save 
write.csv(topmods_tp2_df, file="Outputs/Results/Candidate_models_TP2_heavy_receffort.csv", row.names=FALSE)
avgmod_tp2 <- model.avg(topmods_tp2) 
summary(avgmod_tp2)

## save summary
avg_tp2_summary <- summary(avgmod_tp2, full=F) #pulling out model averages
tp2_summary<-as.data.frame(avg_tp2_summary$coefmat.full) #selecting full model coefficient averages
## get confidence intervals
confint(avgmod_tp2, full=F)
CI <- as.data.frame(confint(avgmod_tp2, full=F)) # get confidence intervals for full model
tp2_summary$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
tp2_summary$CI.max <-CI$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(tp2_summary, keep.rownames = "coefficient") #put rownames into column
names(tp2_summary) <- gsub(" ", "", names(tp2_summary)) # remove spaces from column headers
## get importance values
tp2_summary$parameters <- row.names(tp2_summary)
row.names(tp2_summary) <- 1:nrow(tp2_summary)
importance_tp2 <- MuMIn::sw(avgmod_tp2) ## total precip, temperature and temperature difference in all 34 models
importance_tp2 <- as.data.frame(importance_tp2)
importance_tp2$parameters <- row.names(importance_tp2)
row.names(importance_tp2) <- 1:nrow(importance_tp2)
tp2_summary <- merge(tp2_summary, importance_tp2, by.x="coefficient", by.y="parameters", all=TRUE)
write.csv(tp2_summary, file="Outputs/Results/Average_model_summary_TP2_heavy_receffort.csv", row.names=FALSE)

## Heatmap plots ##
# Interaction between annual temperature & total precipitation
pred.dat_tp2 = expand.grid(temp_diff = median(nmrs_glmm$temp_diff), precip_diff = median(nmrs_glmm$precip_diff),
                           temperature.TP2 = seq(min(nmrs_glmm$temperature.TP2), max(nmrs_glmm$temperature.TP2), length=100),
                           total_precip.TP2 = seq(min(nmrs_glmm$total_precip.TP2), max(nmrs_glmm$total_precip.TP2), length=100),
                           Common_name=unique(nmrs_glmm$Common_name))


# Add predictions
pred.dat_tp2$extinction = predict(avgmod_tp2, newdata=pred.dat_tp2, type="response")
write.csv(pred.dat_tp2, file="Outputs/Results/Predicted_extinction_tp2_heavy.csv", row.names=FALSE)

nmrs_temp_precip <- unique(nmrs_glmm[,c("Hectad","total_precip.TP2", "temperature.TP2")])
precip_temp_TP2 <- ggplot(pred.dat_tp2, aes(total_precip.TP2, temperature.TP2, fill=extinction)) +
  geom_tile() +
  scale_fill_gradient2(low="dodgerblue2", mid="white", high="red",
                       midpoint=0.5) +
  labs(fill="Probability of \nlocal extinction") +
  geom_point(data=nmrs_temp_precip, aes(total_precip.TP2, temperature.TP2), colour="grey", alpha=0.5, inherit.aes = FALSE) +
  xlab(expression("Total precipitation T"[2])) + 
  ylab(expression("Mean temperature T"[2])) + 
  coord_cartesian(xlim = c(quantile(nmrs_glmm$total_precip.TP2, 0.05), quantile(nmrs_glmm$total_precip.TP2, 0.95)),
                  ylim = c(quantile(nmrs_glmm$temperature.TP2, 0.05), quantile(nmrs_glmm$temperature.TP2, 0.95))) +
  theme_classic()
## 1 = persist, 0 = extirpate
## red = high probability of extirpatation
## blue = high probability of persistance
precip_temp_TP2
ggsave(precip_temp_TP2, file="Outputs/Graphs/Temp_precip_TP2_heatmap_heavy.png", height=5, width=7)

## put the 2 interaction plots together with a shared legend

legend <- cowplot::get_legend(precip_temp_TP1 + theme(legend.position = "bottom"))
library(gridExtra)
library(grid)
lay <- rbind(c(1,2), c(3,3))
th <- sum(legend$heights)
final_heatmaps2 <- grid.arrange(precip_temp_TP1 + theme(legend.position="none"), 
                                precip_temp_TP2 + theme(legend.position="none"), legend, 
                                ncol=2, layout_matrix=lay, heights = unit.c(unit(1, "null"), th))
final_heatmaps2
ggsave("Outputs/Graphs/Final_heatmaps_heavy.png", final_heatmaps2, height=4, width=10)
