# Project: Forecasting range shifts of a dioecious plant species under climate change
# Purpose: Map of Common garden experiment for Poa arachnifera across a climatic gradient. 
# Note: Raster files are too large to provide in public repository. They are stored on a local machine
# Authors: Jacob Moutouama, Aldo Compagnoni and Tom Miller
# Date last modified (Y-M-D): 2024-08-02

# remove all objects and clear workspace
rm(list = ls(all=TRUE))

# load packages
library(sp)
library(rgdal)
library(raster)
library(terra)
library(maptools)
library(rgeos)
library(RColorBrewer)
library(tidyverse)
library(dismo)

# Common garden, natural  and source populations locations ---- 
read.csv("https://www.dropbox.com/s/xk4225mn8btqhbm/demography_allsites.csv?dl=1", stringsAsFactors = F) %>% 
  dplyr::select(Latitude,Longitude) %>%
  unique() %>% 
  arrange(Latitude)->garden ## common garden populations

read_rds("https://www.dropbox.com/scl/fi/n3y2mx6xwbnxql84q93v7/poar_ms_quantities.rds?rlkey=fols5060yvrzub8qnqj47nalz&dl=1")$survey_site_table %>% 
  filter(Experimental_source=="yes") %>% 
  dplyr::select(Latitude,Longitude) %>% 
  unique() %>% 
  arrange(Latitude)->source ## source populations

# dir.create("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/occurence")
poar_occ_raw <- gbif(genus="Poa",species="arachnifera",download=TRUE) 
head(poar_occ_raw) 
poar_occ <- subset(poar_occ_raw,(!is.na(lat))&(!is.na(lon))) # here we remove erroneous coordinates, where either the latitude or longitude is missing
cat(nrow(poar_occ_raw)-nrow(poar_occ), "records are removed") # Show the number of records that are removed from the dataset. 
poar_occ %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(Longitude=lon,Latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(Longitude >=-126.374160) &  as.numeric(Longitude <=-95) & country=="United States") %>% 
  unique() %>% 
  arrange(Latitude)->gbif

coordinates(gbif) <- ~ Longitude + Latitude
CRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(gbif) <- CRS1
coordinates(garden) <- ~ Longitude + Latitude
crs(garden) <- CRS1
coordinates(source) <- ~ Longitude + Latitude
crs(source) <- CRS1

# Climatic data past, present and future -----
jacob_path<-"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Future and historical data"  
# tom_path<-"C:/Users/tm9/Dropbox/Miller Lab -- Moutouama/Future and historical data"
choose_path<-jacob_path

## Past -----
raster_1901_1930_list <- list.files(path = paste0(choose_path,"/Historical"),pattern=".tif$",full.names = T)# load bioclimatic layers
clim_1901_1930 <- raster::stack(raster_1901_1930_list) # put all the rasters together
names(clim_1901_1930)<-c("pptdorm","tempdorm","pptgrow","tempgrow") # Change the default names
# plot(clim_1901_1930)
# summary(clim_1901_1930)
Conditions_1901_1930 <-terra::extract(clim_1901_1930,garden)
# head(Conditions_1901_1930)
clim_past<-data.frame(Conditions_1901_1930)

## Present -----
raster_1990_2019_list <- list.files(path = paste0(choose_path,"/Present"),pattern=".tif$",full.names = T)
clim_1990_2019 <- raster::stack(raster_1990_2019_list) # put all rasters together
names(clim_1990_2019)<-c("pptdorm","tempdorm","pptgrow","tempgrow")
Conditions_1990_2019 <- terra::extract(clim_1990_2019,garden)
clim_current<-data.frame(Conditions_1990_2019)
# summary(clim_1990_2019)
# plot(clim_1990_2019)

## Future -----
raster_miroc45_list <- list.files(path = paste0(choose_path,"/Future/MIROC/rcp45"), pattern=".tif$",full.names = T)# load bioclimatic layers
clim_miroc_45 <- terra::rast(raster_miroc45_list) # put all the rasters
clim_miroc_45<-raster::stack(clim_miroc_45)
names(clim_miroc_45)<-c("pptdorm"," tempdorm ","pptgrow","tempgrow") # Change the default names
# plot(clim_miroc_45)
Conditions_miroc85 <- terra::extract(clim_miroc_45,garden)
clim_miroc45<-data.frame(Conditions_miroc85)

raster_miroc85_list <- list.files(path = paste0(choose_path,"/Future/MIROC/rcp85"),pattern=".tif$",full.names = T)# load bioclimatic layers
clim_miroc_85 <- terra::rast(raster_miroc85_list) # put all the raster together
clim_miroc_85<-raster::stack(clim_miroc_85)
names(clim_miroc_85)<-c("pptdorm"," tempdorm ","pptgrow","tempgrow") # Change the defaul
# plot(clim_miroc_85)
Conditions_miroc85 <- terra::extract(clim_miroc_85,garden)
clim_miroc85<-data.frame(Conditions_miroc85)

# Observed data (growing and dormant season temperature and precipitation)----
poar_ppt <- read.csv("https://www.dropbox.com/s/kkga2hf9k1w9ht1/Poa_pr.csv?dl=1", stringsAsFactors = F) # monthly  precipitation 
# head(poar_ppt)
poar_tempmax <- read.csv("https://www.dropbox.com/scl/fi/2eizs4j1sv8g0vgy625qu/tasmax_poar.csv?rlkey=3u6l0765wkmvxulmbk00xuvhk&dl=1", stringsAsFactors = F) # monthly maximum temperature)
poar_tempmin <- read.csv("https://www.dropbox.com/scl/fi/wfrru3bbhxisr6c5yi75y/tasmin_poar.csv?rlkey=y1fihjfwqyq4mfjdkxl4607fs&dl=1", stringsAsFactors = F) # monthly maximum temperature)
poar_temp<-mutate(poar_tempmax,poar_tempmin)
# head(poar_temp)
poar_temp$temp<-(poar_temp$tempmax+poar_temp$tempmin)/2
# head(poar_temp)
poar_ppt$census.year<-ifelse(poar_ppt$month<=5,poar_ppt$year-1,poar_ppt$year) # data were collected in May. Thus, the precipitation for the census year is both the precipitation from June 1 to December of the previous year and precipitation of the current year up to and including May.  
poar_temp$census.year<-ifelse(poar_temp$month<=5,poar_temp$year-1,poar_temp$year)
poar_temp %>% 
  dplyr::select(site,Longitude, Latitude,census.year,temp,month) %>% 
  filter(census.year %in% (2015:2016)) %>%
  mutate(Season=ifelse((month >= 6) & (month <= 9), "dormant", "growing")) %>% 
  group_by(site,census.year,Season) %>%
  dplyr::summarize(site=unique(site),Longitude=unique(Longitude),Latitude=unique(Latitude),temp_season = mean(temp))->poar_tempseason 
# head(poar_tempseason)
poar_ppt %>% 
  dplyr::select(site,Longitude, Latitude,census.year,ppt,month) %>% 
  filter(census.year %in% (2015:2016)) %>%
  mutate(Season=ifelse((month >= 6) & (month <= 9), "dormant", "growing")) %>% 
  group_by(site,census.year,Season) %>%
  dplyr::summarize(site=unique(site),Longitude=unique(Longitude),Latitude=unique(Latitude),ppt_season = sum(ppt))->poar_pptseason 
# head(poar_pptseason)
poar_climseason<- merge(x = poar_tempseason, y = poar_pptseason) # merge temperature and precipitation data
poar_climseason %>% 
  filter(Season=="growing") %>% 
  mutate(tempgrow=temp_season,pptgrow=ppt_season)->poar_growing_2015_2016
poar_climseason %>% 
  filter(Season=="dormant") %>% 
  mutate(tempdorm=temp_season,pptdorm=ppt_season)->poar_dormant_2015_2016
poar_dormant_growing_2015_2016<-mutate(poar_growing_2015_2016,poar_dormant_2015_2016)
poar_dormant_growing_2015_2016 %>% 
  mutate(ztempgrow=scale(tempgrow, scale = TRUE),zpptgrow=scale(pptgrow, scale = TRUE),ztempdorm=scale(tempdorm, scale = TRUE),zpptdorm=scale(pptdorm, scale = TRUE))->poar_2015_2016
# head(poar_2015_2016)
site_cols<-rainbow(length(unique(poar_2015_2016$site)))
year_shapes<-c(1,16)
poar_2015_2016$site_col <- site_cols[as.factor(poar_2015_2016$site)]
poar_2015_2016$year_pch <- year_shapes[as.factor(poar_2015_2016$census.year)]

# Study area shapefile ----
study_area<-terra::vect("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/USA_vector_polygon/States_shapefile.shp")
study_area <- study_area[(study_area$State_Name %in% c("TEXAS","OKLAHOMA","KANSAS")), ]
plot(study_area)

# Clip the climatic rasters
pptdorm_raster<-terra::rast(clim_1990_2019[[1]])
tempdorm_raster<-terra::rast(clim_1990_2019[[2]])
pptgrow_raster<-terra::rast(clim_1990_2019[[3]])
tempgrow_raster<-terra::rast(clim_1990_2019[[4]])
crop_pptdorm_raster <- terra::crop(pptdorm_raster, study_area,mask=TRUE)
crop_tempdorm_raster <- terra::crop(tempdorm_raster, study_area,mask=TRUE)
crop_pptgrow_raster <- terra::crop(pptgrow_raster, study_area,mask=TRUE)
crop_tempgrow_raster <- terra::crop(tempgrow_raster, study_area,mask=TRUE)


# Maps (Figure 1) ----
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/POAR_survey_garden_map.pdf",width=7,height=8,useDingbats = F)
par(mar=c(5,5,1.5,0.5),mfrow=c(2,2))
plot(crop_pptdorm_raster,xlab="Longitude",ylab="Latitude",main="Precipitation (dorm. season)",cex.main=1)
# mtext("A",side = 3, adj = 0,cex=1.25)
plot(study_area,add=T)
plot(gbif,add=T,pch = 23,col="grey50",bg="grey",cex =0.5)
plot(garden,add=T,pch = 21,col="black",cex =1.25,bg="black")
plot(source,add=T,pch = 21,col="black",bg="red",cex =1)
legend(-106.5, 28.5, 
       legend=c( "GBIF occurences","Common garden sites","Source populations"),
       pch = c(23,21,21),
       pt.cex=c(0.55,1.25,1),
       col = c("grey50","black","black"),
       pt.bg=c("grey","black","red"),
       cex = 0.7, 
       bty = "n", 
       horiz = F , 
)

plot(crop_tempdorm_raster,xlab="Longitude",ylab="Latitude",main="Temperature (dorm. season)",cex.main=1)
# mtext("B",side = 3, adj = 0,cex=1.25)
plot(study_area,add=T)
# plot(survey,add=T,pch = 23,col="black",bg="grey",cex =1.5)
plot(gbif,add=T,pch = 23,col="grey50",bg="grey",cex =0.5)
plot(garden,add=T,pch = 21,col="black",cex =1.25,bg="black")
plot(source,add=T,pch = 21,col="black",bg="red",cex =1)

plot(crop_pptgrow_raster,xlab="Longitude",ylab="Latitude",main="Precipitation (grow. season)",cex.main=1)
# mtext("C",side = 3, adj = 0,cex=1.25)
plot(study_area,add=T)
plot(gbif,add=T,pch = 23,col="grey50",bg="grey",cex =0.5)
plot(garden,add=T,pch = 21,col="black",cex =1.25,bg="black")
plot(source,add=T,pch = 21,col="black",bg="red",cex =1)

plot(crop_tempgrow_raster,xlab="Longitude",ylab="Latitude",main="Temperature (grow. season)",cex.main=1)
# mtext("D",side = 3, adj = 0,cex=1.25)
plot(study_area,add=T)
plot(gbif,add=T,pch = 23,col="grey50",bg="grey",cex =0.5)
plot(garden,add=T,pch = 21,col="black",cex =1.25,bg="black")
plot(source,add=T,pch = 21,col="black",bg="red",cex =1)

dev.off()

# Tom's version of the map+climate figure.
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/tom_map_v1.pdf",width=11,height=5)
par(mfrow=c(1,3))
plot(study_area,xlab="Longitude",ylab="Latitude",cex=2,cex.lab=1.2)
plot(gbif,add=T,pch = 23,col="grey50",bg="grey",cex =0.55)
plot(garden,add=T,pch = 21,col="black",cex =2,bg=unique(poar_2015_2016$site_col))
plot(source,add=T,pch = 21,col="black",bg="red",cex =1)
mtext( "A",side = 3, adj = 0,cex=1.25,line=0.2)
legend(-106.5, 28.5, 
       legend=c( "GBIF occurences","Common garden sites","Source populations"),
       pch = c(23,1,21),
       pt.cex=c(0.55,2,1),
       col = c("grey50","black","black"),
       pt.bg=c("grey","white","red"),
       cex = 0.7, 
       bty = "n", 
       horiz = F , 
)

plot(clim_past$pptdorm,clim_past$tempdorm,xlab="Precipitation (mm)",ylab="Temperature (C)",col=unique(poar_2015_2016$site_col),cex=2,xlim=c(150,550),ylim=c(23,34),pch=16,cex.lab=1.2)
mtext( "B",side = 3, adj = 0,cex=1.25)
#points(clim_current$pptdorm,clim_current$tempdorm,col=unique(poar_2015_2016$site_col),cex=1,pch=16)
arrows(clim_past$pptdorm,clim_past$tempdorm,
       clim_current$pptdorm,clim_current$tempdorm,
       length=0.1,col=unique(poar_2015_2016$site_col),lwd=2)
#points(clim_miroc45$pptdorm,clim_miroc45$tempdorm,col=unique(poar_2015_2016$site_col),cex=1.5,pch=1)
#arrows(clim_current$pptdorm,clim_current$tempdorm,
#       clim_miroc45$pptdorm,clim_miroc45$tempdorm,
#       length=0.1,col=unique(poar_2015_2016$site_col),lty=2,lwd=2)
#points(clim_miroc85$pptdorm,clim_miroc85$tempdorm,col=unique(poar_2015_2016$site_col),cex=2,pch=1)
arrows(clim_current$pptdorm,clim_current$tempdorm,
       clim_miroc85$pptdorm,clim_miroc85$tempdorm,
       length=0.1,col=unique(poar_2015_2016$site_col),lwd=2)


plot(clim_past$pptgrow,clim_past$tempgrow,xlab="Precipitation (mm)",ylab="Temperature (C)",col=unique(poar_2015_2016$site_col),cex=2,xlim=c(190,840),ylim=c(7,22),pch=16,cex.lab=1.2)
#points(clim_current$pptgrow,clim_current$tempgrow,col=unique(poar_2015_2016$site_col),cex=1,pch=16)
mtext( "C",side = 3, adj = 0,cex=1.25)
arrows(clim_past$pptgrow,clim_past$tempgrow,
       clim_current$pptgrow,clim_current$tempgrow,
       length=0.1,col=unique(poar_2015_2016$site_col),lwd=2)
#points(clim_miroc45$pptgrow,clim_miroc45$tempgrow,col=unique(poar_2015_2016$site_col),cex=1.5,pch=1)
#arrows(clim_current$pptgrow,clim_current$tempgrow,
#       clim_miroc45$pptgrow,clim_miroc45$tempgrow,
#       length=0.1,col=unique(poar_2015_2016$site_col),lty=2,lwd=2)
#points(clim_miroc85$pptgrow,clim_miroc85$tempgrow,col=unique(poar_2015_2016$site_col),cex=2,pch=1)
arrows(clim_current$pptgrow,clim_current$tempgrow,
       clim_miroc85$pptgrow,clim_miroc85$tempgrow,
       length=0.1,col=unique(poar_2015_2016$site_col),lwd=2)
dev.off()


pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/tom_map_v2.pdf",width=11,height=5)

par(mfrow=c(1,3))
plot(study_area,xlab="Longitude",ylab="Latitude",cex.lab=1.2)
plot(gbif,add=T,pch = 23,col="grey50",bg="grey",cex =0.55)
plot(garden,add=T,pch = 3,col="black",cex =3,bg=unique(poar_2015_2016$site_col))
plot(source,add=T,pch = 21,col="black",bg="red",cex =1)
mtext( "A",side = 3, adj = 0,cex=1.25,line=0.2)
legend(-106, 28.5, 
       legend=c( "GBIF occurences","Common garden sites","Source populations"),
       pch = c(23,3,21),
       pt.cex=c(0.55,2,1),
       col = c("grey50","black","black"),
       pt.bg=c("grey","black","red"),
       cex = 0.9, 
       bty = "n", 
       horiz = F , 
)

plot(clim_past$pptdorm,clim_past$tempdorm,xlab="Precipitation (mm)",ylab="Temperature (C)",col=alpha("black",0.25),cex=2,xlim=c(150,550),ylim=c(23,34),pch=16,cex.lab=1.2)
points(clim_miroc45$pptdorm,clim_miroc45$tempdorm,col=alpha("blue",0.25),cex=2,pch=16)
points(clim_miroc85$pptdorm,clim_miroc85$tempdorm,col=alpha("red",0.25),cex=2,pch=16)
points(poar_2015_2016$pptdorm,poar_2015_2016$tempdorm,col="black",cex=2,pch=3)
mtext( "B",side = 3, adj = 0,cex=1.25)
title(main="Dormant season",adj=0.5,line=0.5)
legend("topright",legend=c("Past (1901-1930)","Observed (2014-2016)","RCP4.5 (2071-2100)","RCP8.5 (2071-2100)"),pch=c(16,3,16,16),col=c(alpha("black",0.25),"black",alpha("blue",0.25),alpha("red",0.25)),cex = 0.9)
#points(clim_current$pptdorm,clim_current$tempdorm,col="black",cex=1,pch=16)


plot(clim_past$pptgrow,clim_past$tempgrow,xlab="Precipitation (mm)",ylab="Temperature (C)",col=alpha("black",0.25),cex=2,xlim=c(190,840),ylim=c(7,22),pch=16,cex.lab=1.2)
points(clim_miroc45$pptgrow,clim_miroc45$tempgrow,col=alpha("blue",0.25),cex=2,pch=16)
points(clim_miroc85$pptgrow,clim_miroc85$tempgrow,col=alpha("red",0.25),cex=2,pch=16)
points(poar_2015_2016$pptgrow,poar_2015_2016$tempgrow,col="black",cex=2,pch=3)
mtext( "C",side = 3, adj = 0,cex=1.25)
title(main="Growing season",adj=0.5,line=0.5)
#points(clim_current$pptgrow,clim_current$tempgrow,col="black",cex=1,pch=16)
dev.off()

