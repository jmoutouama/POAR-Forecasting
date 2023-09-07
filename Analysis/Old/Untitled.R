# Code to dowload US corrdianate to project the model on the geograohical space. 
library(sf)
library(ggmap)
outline_map <- map_data("world")
states_shape <- map_data("state")
counties_shape <- map_data("county")
US_coordinate<-data.frame(long=states_shape$long,lat=states_shape$lat,region=states_shape$region)
write.csv(US_coordinate,"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/US_coordinate.csv",row.names = F)
