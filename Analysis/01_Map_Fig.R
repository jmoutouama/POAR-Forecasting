#### Project: Thunbergia atacorensis demography 2019-2021
#### PURPOSE: Create map of T. atacorensis study populations (Fig. 1a)
#### Note: GIS files are too large to provide in public repository
#### AUTHOR: Jacob Moutouama
#### DATE LAST MODIFIED: 20210909

# remove objects and clear workspace
rm(list = ls(all=TRUE))

#load packages

library(tmap)
library(sf)
library(tidyverse)
library(geojsonio)
library(sp)
library(rgeos)
library(ggmap)
library(maptools)



#*************************************
# Plot map of demography populations, Fig. 2a
#*************************************

#read in T atacorensis demography data and extract lat/lon for each population
Thunbergia <- read.csv("/Users/jmoutouama/Dropbox/PhD Project/Chapter3/Demography Thunbergia/Data/occurenceplot.csv", header = T)

#read the shapefile of the study area
states_shape <- map_data("state")
US = states_shape %>% filter(!region %in% c("Alaska", "Hawaii", "Puerto Rico"))


#project demography site locations

Thunbergia_sf<-st_as_sf(Thunbergia,coords = c("longitude","latitude"), crs = 4326,agr = "identity")

#plot the map with density per plot
tm_shape(District) +
  tm_polygons(col = "grey",border.col = "grey")+
  tm_shape(Thunbergia_sf)+
  tm_dots(size="Density.per.plot",col = "black", title.size = "Density/plot")+
  tmap::tm_compass(type = "4star", position = c("left", "top")) +
  tmap::tm_scale_bar(breaks = c(0, 45, 90), text.size = 0.75,position = c("left", "bottom"))+
  tmap::tm_graticules(lines=FALSE,labels.rot=c(0,90))

#plot the map without density per plot
tm_shape(District) +
  tm_polygons(col = "grey",border.col = "grey50")+
  tm_shape(Thunbergia_sf)+
  tm_dots(size=0.2,col = "black", title.size = "Plot")+
  tmap::tm_compass(type = "8star", position = c("right", "top")) +
  tmap::tm_scale_bar(breaks = c(0, 45, 90), text.size = 0.5,position = c("left", "bottom"))+
  tmap::tm_graticules(lines=FALSE,labels.rot=c(0,90))



# 

tm_shape(Benin) +
  tm_polygons(col = "grey",border.col = "grey50")






