# PROJECT: Forecasting range shifts of a dioecious plant species under climate change
# PURPOSE: Map of Common garden experiment for Poa arachnifera across a climatic gradient. 
## Note: Raster files are too large to provide in public repository. They are stored on a local machine
### AUTHOR: Jacob Moutouama, Aldo Compagnoni and Tom Miller
#### DATE LAST MODIFIED (Y-M-D): 2024-08-02

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

# Climatic data (1990-2019) -----
jacob_path<-"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Future and historical data"  
# tom_path<-"C:/Users/tm9/Dropbox/Miller Lab -- Moutouama/Future and historical data"
choose_path<-jacob_path
raster_1990_2019_list <- list.files(path = paste0(choose_path,"/Present"),pattern=".tif$",full.names = T)
clim_1990_2019 <- raster::stack(raster_1990_2019_list) # put all rasters together
names(clim_1990_2019)<-c("pptdorm","tempdorm","pptgrow","tempgrow")
# summary(clim_1990_2019)
# plot(clim_1990_2019)

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





