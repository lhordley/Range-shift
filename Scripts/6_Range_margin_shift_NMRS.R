##########################
#### user: Lisbeth Hordley
#### date: July 2021
#### info: Range margin shift analysis using NMRS data

rm(list = ls())
options(scipen=999)

## load packages
library(data.table)
library(tidyverse)
library(lubridate)
library(sf)
library(ggplot2)
library(plyr)
library(ggpubr)
library(Rmisc)

## NMRS data for cool-adapted moths with elevation and temperature data at a 10km scale
nmrsdata <- readRDS("Data/NMRS/NMRS_cool_moths_final.rds") ## NMRS data for all hectads and all years with elevation
## hectad recording levels
hec_records <- read.csv("Data/NMRS/Hectad_recording_levels_1975_1991_2012_2016.csv", header=TRUE)

###############################################################
###################### RECORDED HECTADS #######################
###############################################################

## filter to recorded hectads only (i.e. those recorded once in both time periods)
hectads <- unique(hec_records$HECTAD) ## 2191
nmrsdata_rec <- nmrsdata[which(nmrsdata$Hectad %in% hectads), ]
length(unique(nmrsdata_rec$Hectad)) ## 1886 hectads - fewer because the cool-adapted species aren't found in all hectads
length(unique(nmrsdata_rec$Common_name)) ## 77 species

### Filter: species need to occupy at least 20 hectads in each time period to calculate range shifts
# count the number of hectads per species for each year
nmrsdata_rec$PRESENCE <- 1
tp_species_hecs <- ddply(nmrsdata_rec, .(Scientific_name,Common_name, Time_period), summarise,
                         HECTADS = sum(PRESENCE)) ## 154 (77 species in two time periods)
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_hecs <- tp_species_hecs %>% group_by(Common_name) %>% filter(all(HECTADS>=20))
species_keep <- unique(tp_species_hecs$Common_name) ## 51 species (lose 26)
## mostly removes rare/restricted species
nmrsdata_rec <-  nmrsdata_rec[which(nmrsdata_rec$Common_name %in% species_keep),] 
length(unique(nmrsdata_rec$Hectad)) ## 1883 hectads
length(unique(nmrsdata_rec$Common_name)) ## 51 species

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
low_10_elev_rec_hecs <- nmrsdata_rec %>% group_by(Common_name, Time_period) %>% arrange(elevation10x10km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_rec_hecs <- summarySE(low_10_elev_rec_hecs, measurevar="elevation10x10km", groupvars=c("Common_name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_rec_hecs, 
          elevation10x10km[Time_period == "1975-1991"] - elevation10x10km[Time_period == "2012-2016"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = <0.001 - distribution is NOT normal
qqnorm(d)
qqline(d) ## not good

## wilcoxon signed-rank test instead
low_elev_rec <- wilcox.test(elevation10x10km ~ Time_period, data = low_elev_rec_hecs, paired = TRUE)
low_elev_rec ## significant (p=0.007) - elevation is higher in TP1 compared to TP2

## plot result
rec_low_elev <- ggpaired(low_elev_rec_hecs, x = "Time_period", y = "elevation10x10km",
                         color = "Time_period", line.color = "gray", line.size = 0.4,
                         palette = "jco", id="Common_name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
rec_low_elev
## low elevation boundary moving downhill over time
ggsave(rec_low_elev, file="Graphs/Low_elevation_shift_rec_hecs.png")

################# LOW LATITUDE BOUNDARY SHIFT - mean latitude of 10 most southerly hectads for each time period

## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_rec_hecs <- nmrsdata_rec %>% group_by(Common_name, Time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_rec_hecs <- summarySE(low_10_lat_rec_hecs, measurevar="lat", groupvars=c("Common_name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_lat_rec_hecs, 
          lat[Time_period == "1975-1991"] - lat[Time_period == "2012-2016"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value <0.001 - non-normal distribution
qqnorm(d)
qqline(d) ## looks ok

## wilcoxon signed-rank test instead
low_lat_rec <- wilcox.test(lat ~ Time_period, data = low_lat_rec_hecs, paired = TRUE)
low_lat_rec ## significant

## plot anyway
rec_low_lat <- ggpaired(low_lat_rec_hecs, x = "Time_period", y = "lat",
                        color = "Time_period", line.color = "gray", line.size = 0.4,
                        palette = "jco", id="Common_name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
rec_low_lat
ggsave(rec_low_lat, file="Graphs/Low_latitude_shift_rec_hecs.png")
##low latitude boundary moving northwards
## Sword grass shows massive northwards shift



##################################################################
##################### WELL RECORDED HECTADS ######################
##################################################################

## remove Lunar thorn
nmrsdata <- nmrsdata[!nmrsdata$Common_name=="Lunar Thorn",]

## filter to recorded hectads only (i.e. those recorded once in both time periods)
well_heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded" | RECORDING.LEVEL=="Well recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
well_heavy_hecs <- well_heavy_hecs[which(well_heavy_hecs$n_row > 1), ]
well_heavy_hecs$n_row <- NULL ## 1084 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_well <- nmrsdata[which(nmrsdata$Hectad %in% well_heavy_hecs$HECTAD), ]
length(unique(nmrsdata_well$Hectad)) ## 982 hectads (639 before - gained 343 hectads)
length(unique(nmrsdata_well$Common_name)) ## 77 species - don't lose any species!

########## Range margin shifts
####### Filter: species need to occupy at least 20 hectads in each time period to calculate range shifts
# count the number of hectads per species for each year
nmrsdata_well$PRESENCE <- 1
tp_species_hecs <- ddply(nmrsdata_well, .(Scientific_name,Common_name, Time_period), summarise,
                         HECTADS = sum(PRESENCE)) ## 144 (72 species in two time periods)
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_hecs <- tp_species_hecs %>% group_by(Common_name) %>% filter(all(HECTADS>=20))
species_keep <- unique(tp_species_hecs$Common_name) ## 46 species
## mostly removes rare/restricted species
nmrsdata_well <-  nmrsdata_well[which(nmrsdata_well$Common_name %in% species_keep),] 
## remove duplicated records at 10km resolution - same species/year/hectad (probably have different 1km sites)
length(unique(nmrsdata_well$Hectad)) ## 982 hectads
length(unique(nmrsdata_well$Common_name)) ## 46 species

## map
worldmap = map_data('world')
well_heavy_hecs <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  geom_point(data =nmrsdata_well, 
             aes(x = lon, y=lat), size=1) +
  theme(title = element_text(size = 12))
well_heavy_hecs
## save graph
#ggsave(well_heavy_hecs, file="../../Maps/Upland_moth_16spp_well_heavy_hecs.png")

################# LOW ELEVATION BOUNDARY SHIFT - mean elevation of 10 lowest elevation hectads for each time period

## order elevation by ascending value
## take the top 10 values 
low_10_elev_wh_hecs <- nmrsdata_well %>% group_by(Common_name, Time_period) %>% arrange(elevation10x10km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_wh_hecs <- summarySE(low_10_elev_wh_hecs, measurevar="elevation10x10km", groupvars=c("Common_name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_wh_hecs, 
          elevation10x10km[Time_period == "1975-1991"] - elevation10x10km[Time_period == "2012-2016"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = <0.001 - distribution is NOT normal
qqnorm(d)
qqline(d) ## not good

## wilcoxon signed-rank test instead
low_elev_wh <- wilcox.test(elevation10x10km ~ Time_period, data = low_elev_wh_hecs, paired = TRUE)
low_elev_wh ## non-significant

## plot anyway
low_elev_wh <- ggpaired(low_elev_wh_hecs, x = "Time_period", y = "elevation10x10km",
                        color = "Time_period", line.color = "gray", line.size = 0.4,
                        palette = "jco", id="Common_name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
low_elev_wh
ggsave(low_elev_wh, file="Graphs/Low_elevation_shift_well_hecs.png")

plot <- ggplot(low_elev_wh_hecs, aes(x=Time_period, y=elevation10x10km)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = "black"))
plot


################# LOW LATITUDE BOUNDARY SHIFT - mean latitude of 10 most southerly hectads for each time period

## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_wh_hecs <- nmrsdata_well %>% group_by(Common_name, Time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_wh_hecs <- summarySE(low_10_lat_wh_hecs, measurevar="lat", groupvars=c("Common_name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_lat_wh_hecs, 
          lat[Time_period == "1975-1991"] - lat[Time_period == "2012-2016"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value <0.001 - non-normal distribution
qqnorm(d)
qqline(d) ## not good

## wilcox.test
low_lat_wh <- wilcox.test(lat ~ Time_period, data = low_lat_wh_hecs, paired = TRUE, exact=F, conf.int = T, conf.level = .95)
low_lat_wh ## significant!

df2 <- low_lat_wh_hecs %>% 
  pivot_wider(id_cols = Common_name, 
              values_from = lat, 
              names_from = Time_period)
df2$diff <- df2$`2012-2016` - df2$`1975-1991`
df2$cat <- ifelse(df2$diff>= 0, "north", "south")
table(df2$cat)

library(showtext)
font_add_google("Montserrat", "Montserrat")
font_add_google("Roboto", "Roboto")
font_files()    
font_add("Calibri", "calibri.ttf")
myFont1 <- "Calibri"

install.packages("showtext")
library(showtext)
list_of_fonts <- as.data.frame(font_files())

grep(pattern = "Calibri", x = list_of_fonts$family, ignore.case = TRUE, value = TRUE)

worldmap = map_data('world')
worldmap <- worldmap[!worldmap$region=="Ireland",]
worldmap <- worldmap[!worldmap$subregion=="Northern Ireland",]
poster_margin <- ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  theme_void() + 
  #geom_point(data =nmrsdata_well[nmrsdata_well$Common_name=="Small Autumnal Moth",], 
             #aes(x = lon, y=lat, colour=Time_period), size=1.7, alpha=0.3) +
  #scale_color_manual(values = c("1975-1991" = "#D55E00", "2012-2016" = "#0072B2")) +
  geom_segment(aes(x=-6,xend=2.3,y=51.94290,yend=51.94290), colour="#D55E00", lwd=2, linetype=5) +
  geom_segment(aes(x=-6,xend=2.3,y=53.32317,yend=53.32317), colour="#0072B2", lwd=2, linetype=5) +
  theme(legend.position="none", plot.margin=unit(c(0,0,0,0), "null")) 
poster_margin
ggsave(poster_margin, file="Graphs/poster_method1.png", height=20, width=14, units="cm")

## plot anyway
low_lat_wh <- ggpaired(low_lat_wh_hecs, x = "Time_period", y = "lat",
                       color = "Time_period", line.color = "gray", line.size = 0.4,
                       palette = "jco", id="Common_name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="wilcox", paired = TRUE)
low_lat_wh
ggsave(low_lat_wh, file="Graphs/Low_latitude_shift_well_hecs.png")
## significant northwards shift in cool-adapted species

plot2 <- ggplot(low_lat_wh_hecs, aes(x=Time_period, y=lat, colour=Time_period)) +
  scale_colour_manual(values=c("#D55E00", "#0072B2")) +
  labs(x="", y="") +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = "black"))
plot2



#####################################################################
##################### HEAVILY RECORDED HECTADS ######################
#####################################################################

## remove Lunar thorn
nmrsdata <- nmrsdata[!nmrsdata$Common_name=="Lunar Thorn",]

## filter to recorded hectads only (i.e. those recorded once in both time periods)
heavy_hecs <- hec_records %>%
  filter(RECORDING.LEVEL == "Heavily recorded") %>%
  group_by(HECTAD) %>%
  dplyr::summarise(n_row=n())
heavy_hecs <- heavy_hecs[which(heavy_hecs$n_row > 1), ]
heavy_hecs$n_row <- NULL ## 1084 hectads
## merge into nmrs data and filter by years so only have years within time periods
nmrsdata_heavy <- nmrsdata[which(nmrsdata$Hectad %in% heavy_hecs$HECTAD), ]
length(unique(nmrsdata_heavy$Hectad)) ## 592
length(unique(nmrsdata_heavy$Common_name)) ## 72 species

########## Range margin shifts
####### Filter: species need to occupy at least 20 hectads in each time period to calculate range shifts
# count the number of hectads per species for each year
nmrsdata_heavy$PRESENCE <- 1
tp_species_hecs <- ddply(nmrsdata_heavy, .(Scientific_name,Common_name, Time_period), summarise,
                         HECTADS = sum(PRESENCE)) ## 144 (72 species in two time periods)
## only keep species if they occur in at least 20 hectads in BOTH time periods
tp_species_hecs <- tp_species_hecs %>% group_by(Common_name) %>% filter(all(HECTADS>=20))
species_keep <- unique(tp_species_hecs$Common_name) ## 46 species
## mostly removes rare/restricted species
nmrsdata_heavy <-  nmrsdata_heavy[which(nmrsdata_heavy$Common_name %in% species_keep),] 
## remove duplicated records at 10km resolution - same species/year/hectad (probably have different 1km sites)
length(unique(nmrsdata_heavy$Hectad)) ## 591 hectads
length(unique(nmrsdata_heavy$Common_name)) ## 38 species

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
#ggsave(well_heavy_hecs, file="../../Maps/Upland_moth_16spp_well_heavy_hecs.png")

################# LOW ELEVATION BOUNDARY SHIFT - mean elevation of 10 lowest elevation hectads for each time period

## order elevation by ascending value
## take the top 10 values 
low_10_elev_heavy_hecs <- nmrsdata_heavy %>% group_by(Common_name, Time_period) %>% arrange(elevation10x10km) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_elev_heavy_hecs <- summarySE(low_10_elev_heavy_hecs, measurevar="elevation10x10km", groupvars=c("Common_name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_elev_heavy_hecs, 
          elevation10x10km[Time_period == "1975-1991"] - elevation10x10km[Time_period == "2012-2016"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = <0.001 - distribution is NOT normal
qqnorm(d)
qqline(d) ## not good

## wilcoxon signed-rank test instead
low_elev_wh <- wilcox.test(elevation10x10km ~ Time_period, data = low_elev_heavy_hecs, paired = TRUE)
low_elev_wh ## significant

## plot anyway
low_elev_wh <- ggpaired(low_elev_heavy_hecs, x = "Time_period", y = "elevation10x10km",
                        color = "Time_period", line.color = "gray", line.size = 0.4,
                        palette = "jco", id="Common_name")+
  xlab("Time period")+
  ylab("Low elevation boundary")+
  stat_compare_means(method="wilcox.test", paired = TRUE)
low_elev_wh
ggsave(low_elev_wh, file="Graphs/Low_elevation_shift_heavy_hecs.png")
## significant shift uphill over time
plot <- ggplot(low_elev_heavy_hecs, aes(x=Time_period, y=elevation10x10km)) +
  geom_boxplot() +
  theme_classic()
plot


################# LOW LATITUDE BOUNDARY SHIFT - mean latitude of 10 most southerly hectads for each time period

## order latitude by ascending value
## take the top 10 values & average those

low_10_lat_heavy_hecs <- nmrsdata_heavy %>% group_by(Common_name, Time_period) %>% arrange(lat) %>% slice(seq_len(10))
## now calculate mean of the 10 lowest hectads + summary stats
low_lat_heavy_hecs <- summarySE(low_10_lat_heavy_hecs, measurevar="lat", groupvars=c("Common_name","Time_period"))

## check for normality 
# compute the difference
d <- with(low_lat_heavy_hecs, 
          lat[Time_period == "1975-1991"] - lat[Time_period == "2012-2016"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value <0.001 - non-normal distribution
qqnorm(d)
qqline(d) ## not good

## wilcox.test
low_lat_wh <- wilcox.test(lat ~ Time_period, data = low_lat_heavy_hecs, paired = TRUE)
low_lat_wh ## significant

## plot anyway
low_lat_wh <- ggpaired(low_lat_heavy_hecs, x = "Time_period", y = "lat",
                       color = "Time_period", line.color = "gray", line.size = 0.4,
                       palette = "jco", id="Common_name")+
  xlab("Time period")+
  ylab("Low latitude boundary")+
  stat_compare_means(method="wilcox", paired = TRUE)
low_lat_wh
ggsave(low_lat_wh, file="Graphs/Low_latitude_shift_heavy_hecs.png")
## significant northwards shift in cool-adapted species

plot2 <- ggplot(low_lat_heavy_hecs, aes(x=Time_period, y=lat)) +
  geom_boxplot() +
  theme_classic()
plot2
