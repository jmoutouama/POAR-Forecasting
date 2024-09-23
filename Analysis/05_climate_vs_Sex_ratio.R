# Project: Forecasting range shifts of a dioecious plant species under climate change
# Question: What is the impact of climate change on operational sex ratio throughout the range? 
# Note: Raster files are too large to provide in public repository. They are stored on a local machine
# Authors: Jacob Moutouama, Aldo Compagnoni and Tom Miller
# Date last modified (Y-M-D): 2024-08-03

# remove all objects and clear workspace
rm(list = ls(all=TRUE))
# load packages
library(rstan)
# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(13)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
options(tidyverse.quiet = TRUE)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(bayesplot)
# install.packages("countreg",repos = "http://R-Forge.R-project.org")
#library(countreg)
library(rmutil)
library(actuar)
#library(SPEI)
library(LaplacesDemon)
library(ggpubr)
library(raster)
library(rgdal)
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("scater")
# library(scater)
library(BiocManager)
library(swfscMisc)
library(do)
library(tidyterra)

# Define some basic functions that we'll use later 
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply(deparse )
}

invlogit<-function(x){exp(x)/(1+exp(x))}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Climatic data past, present and future -----
jacob_path<-"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Future and historical data"  
# tom_path<-"C:/Users/tm9/Dropbox/Miller Lab -- Moutouama/Future and historical data"
choose_path<-jacob_path

## Common garden populations locations  
read.csv("https://www.dropbox.com/s/xk4225mn8btqhbm/demography_allsites.csv?dl=1", stringsAsFactors = F) %>% 
  dplyr::select(Latitude,Longitude) %>%
  unique() %>% 
  arrange(Latitude)->garden ## common garden populations

coordinates(garden) <- ~ Longitude + Latitude
CRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(garden) <- CRS1 

## Past -----
raster_1901_1930_list <- list.files(path = paste0(choose_path,"/Historical"),pattern=".tif$",full.names = T)# load bioclimatic layers
clim_1901_1930 <- raster::stack(raster_1901_1930_list) # put all the rasters together
names(clim_1901_1930)<-c("pptdorm","tempdorm","pptgrow","tempgrow") # Change the default names
plot(clim_1901_1930)
Conditions_1901_1930 <-terra::extract(clim_1901_1930,garden)
clim_past<-data.frame(Conditions_1901_1930)

## Present -----
raster_1990_2019_list <- list.files(path = paste0(choose_path,"/Present"),pattern=".tif$",full.names = T)
clim_1990_2019 <- raster::stack(raster_1990_2019_list) # put all rasters together
names(clim_1990_2019)<-c("pptdorm","tempdorm","pptgrow","tempgrow")
Conditions_1990_2019 <- terra::extract(clim_1990_2019,garden)
clim_current<-data.frame(Conditions_1990_2019)

## Future -----
### MIROC5----
raster_miroc45_list <- list.files(path = paste0(choose_path,"/Future/MIROC/rcp45"), pattern=".tif$",full.names = T)# load bioclimatic layers
clim_miroc_45 <- terra::rast(raster_miroc45_list) # put all the rasters
clim_miroc_45<-raster::stack(clim_miroc_45)
names(clim_miroc_45)<-c("pptdorm"," tempdorm ","pptgrow","tempgrow") # Change the default names
Conditions_miroc85 <- terra::extract(clim_miroc_45,garden)
clim_miroc45<-data.frame(Conditions_miroc85)
raster_miroc85_list <- list.files(path = paste0(choose_path,"/Future/MIROC/rcp85"),pattern=".tif$",full.names = T)# load bioclimatic layers
clim_miroc_85 <- terra::rast(raster_miroc85_list) # put all the raster together
clim_miroc_85<-raster::stack(clim_miroc_85)
names(clim_miroc_85)<-c("pptdorm"," tempdorm ","pptgrow","tempgrow") # Change the defaul
Conditions_miroc85 <- terra::extract(clim_miroc_85,garden)
clim_miroc85<-data.frame(Conditions_miroc85)

### ACCESS----
raster_access_list_45 <- list.files(path = paste0(choose_path,"/Future/ACCESS/rcp45"),pattern=".tif$",full.names = T)# load bioclimatic layers
clim_access_45<-raster::stack(raster_access_list_45)
names(clim_access_45)<-c("pptdorm"," tempdorm ","pptgrow","tempgrow") # Change the default names
plot(clim_access_45)
Conditions_access_45 <- terra::extract(clim_access_45,garden)
clim_access45<-data.frame(Conditions_access_45)
raster_access_list_85 <- list.files(path = paste0(choose_path,"/Future/ACCESS/rcp85"),pattern=".tif$",full.names = T)# load bioclimatic layers
clim_access_85<-raster::stack(raster_access_list_85)
names(clim_access_85)<-c("pptdorm"," tempdorm ","pptgrow","tempgrow") # Change the default names
plot(clim_access_85)
Conditions_access_85 <- terra::extract(clim_access_85,garden)
clim_access85<-data.frame(Conditions_access_85)

### CMCC ----
raster_cmcc_list_45 <- list.files(path = paste0(choose_path,"/Future/CMCC/rcp45"),pattern=".tif$",full.names = T)# load bioclimatic layers
clim_cmcc_45<-raster::stack(raster_cmcc_list_45)
names(clim_cmcc_45)<-c("pptdorm"," tempdorm ","pptgrow","tempgrow") # Change the default names
plot(clim_cmcc_45)
Conditions_cmcc_45 <- terra::extract(clim_cmcc_45,garden)
clim_cmcc45<-data.frame(Conditions_cmcc_45)
raster_cmcc_list_85 <- list.files(path = paste0(choose_path,"/Future/CMCC/rcp85"),pattern=".tif$",full.names = T)# load bioclimatic layers
clim_cmcc_85<-raster::stack(raster_cmcc_list_85)
names(clim_cmcc_85)<-c("pptdorm"," tempdorm ","pptgrow","tempgrow") # Change the default names
plot(clim_cmcc_85)
Conditions_cmcc_85 <- terra::extract(clim_cmcc_85,garden)
clim_cmcc85<-data.frame(Conditions_cmcc_85)

### CESM1----
raster_cesm_list_45 <- list.files(path = paste0(choose_path,"/Future/CESM1/rcp45"),pattern=".tif$",full.names = T)# load bioclimatic layers
clim_cesm_45<-raster::stack(raster_cesm_list_45)
names(clim_cesm_45)<-c("pptdorm"," tempdorm ","pptgrow","tempgrow") # Change the default names
plot(clim_cesm_45)
Conditions_cesm_45 <- terra::extract(clim_cesm_45,garden)
clim_cesm45<-data.frame(Conditions_cesm_45)
raster_cesm_list_85 <- list.files(path = paste0(choose_path,"/Future/CESM1/rcp85"),pattern=".tif$",full.names = T)# load bioclimatic layers
clim_cesm_85<-raster::stack(raster_cesm_list_85)
names(clim_cesm_85)<-c("pptdorm"," tempdorm ","pptgrow","tempgrow") # Change the default names
plot(clim_cesm_85)
Conditions_cesm_85 <- terra::extract(clim_cesm_85,garden)
clim_cesm85<-data.frame(Conditions_cesm_85)

## Observed data (growing and dormant season temperature and precipitation) ----
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

# Demography data ----
poar_allsites <- read.csv("https://www.dropbox.com/s/xk4225mn8btqhbm/demography_allsites.csv?dl=1", stringsAsFactors = F)#common garden data
poar_allsites$census.year<-poar_allsites$year-1 # Add census year to match with climate data
# unique(poar_allsites$census.year)
poar_allsites %>% 
  dplyr::select(everything()) %>% 
  filter(census.year %in% (2015:2016))-> poar_allsites_2015_2016 # Drop the first census to match with the seasonal model (growing and dormant season temperature and precipitation)

# Matching seasonal data with observed demographic data
poar.clim_seasonal <- left_join(x = poar_allsites_2015_2016,y =poar_2015_2016,by=c("site","census.year","Longitude","Latitude")) # merge the demographic data with climatic data for each site
# head(poar.clim_seasonal)

## Read and format data to build the model
# Survival data
poar.clim_seasonal %>% 
  subset( tillerN_t0 > 0 )%>%
  dplyr::select(census.year, Code, site, Block, Sex, 
                Longitude, Latitude, 
                tillerN_t0, surv_t1,zpptgrow,zpptdorm,ztempgrow,ztempdorm,site)%>% 
  na.omit %>% 
  mutate( site= site %>% as.factor %>% as.numeric,
          Block = Block %>% as.factor %>% as.numeric,
          Sex = Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric ) %>%
  mutate( log_size_t0   = log(tillerN_t0),
          z_ppt_t1_grow = zpptgrow,
          z_ppt_t1_dorm = zpptdorm,
          z_temp_t1_grow = ztempgrow,
          z_temp_t1_dorm = ztempdorm)->poarclim_seasonal.surv
# growth data
poar.clim_seasonal %>% 
  subset( tillerN_t0 > 0 & tillerN_t1 > 0)%>%
  dplyr::select( census.year, Code, site, Block, Sex, 
                 Longitude, Latitude, 
                 tillerN_t0, tillerN_t1,surv_t1,zpptgrow,zpptdorm,ztempgrow,ztempdorm,site)%>% 
  na.omit %>% 
  mutate( site= site %>% as.factor %>% as.numeric,
          Block = Block %>% as.factor %>% as.numeric,
          Sex= Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric ) %>%
  mutate( log_size_t0   = log(tillerN_t0),
          z_ppt_t1_grow = zpptgrow,
          z_ppt_t1_dorm = zpptdorm,
          z_temp_t1_grow = ztempgrow,
          z_temp_t1_dorm = ztempdorm)->poarclim_seasonal.grow
# flowering data
poar.clim_seasonal %>% 
  subset( tillerN_t1 > 0 )%>%
  dplyr::select(census.year, Code, site, Block, Sex, 
                Longitude, Latitude, 
                tillerN_t1, flowerN_t1,flow_t1,zpptgrow,zpptdorm,ztempgrow,ztempdorm,site)%>% 
  na.omit %>% 
  mutate( site= site %>% as.factor %>% as.numeric,
          Block = Block %>% as.factor %>% as.numeric,
          Sex= Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric ) %>%
  mutate( log_size_t1   = log(tillerN_t1),
          z_ppt_t1_grow = zpptgrow ,
          z_ppt_t1_dorm = zpptdorm ,
          z_temp_t1_grow = ztempgrow ,
          z_temp_t1_dorm = ztempdorm)->poarclim_seasonal.flow
# panicle data
poar.clim_seasonal %>% 
  subset( flowerN_t1 > 0 & tillerN_t1 > 0 )%>%
  dplyr::select( census.year, Code, site, Block, Sex, 
                 Longitude, Latitude, 
                 tillerN_t1, flowerN_t1,zpptgrow,zpptdorm,ztempgrow,ztempdorm,site)%>% 
  na.omit %>% 
  mutate( site = site %>% as.factor %>% as.numeric,
          Block = Block %>% as.factor %>% as.numeric,
          Sex= Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric,
          panic_t1 = flowerN_t1) %>%
  mutate( 
    log_size_t1 = log(tillerN_t1),
    z_ppt_t1_grow = zpptgrow,
    z_ppt_t1_dorm = zpptdorm ,
    z_temp_t1_grow = ztempgrow ,
    z_temp_t1_dorm = ztempdorm )->poarclim_seasonal.panic
#seed viability and germination
viabVr <- read.csv("https://www.dropbox.com/s/jfkgoxgv8o1fgqx/viability.csv?dl=1")
# seed viability
viabVr %>% 
  dplyr::select( plot, totS, yesMaybe, sr_f ) %>% 
  dplyr::rename( SR        = sr_f,
                 y_viab = yesMaybe,
                 tot_seeds_viab = totS) %>% 
  dplyr::select(y_viab, tot_seeds_viab, SR ) %>% 
  na.omit ->viab
# seed germination
viabVr %>% 
  dplyr::select( plot, germTot, germFail, sr_f ) %>% 
  dplyr::rename( SR        = sr_f,
                 y_germ    = germTot ) %>% 
  mutate(tot_seeds_germ = y_germ + germFail ) %>% 
  dplyr::select(y_germ, tot_seeds_germ, SR ) %>% 
  na.omit->germ
# seeds per panicle
viabVr %>% 
  dplyr::select(SeedN)  %>% 
  na.omit->seeds


# Sex ratio for observed data----
garden_osr <- poar.clim_seasonal %>% 
  dplyr::select(year, site, Sex,zpptgrow,zpptdorm,ztempgrow,ztempdorm,flowerN_t1) %>% 
  group_by(year,zpptgrow,zpptdorm,ztempgrow,ztempdorm,site, Sex) %>% 
  summarise(total_panicles = sum(flowerN_t1,na.rm=T)) %>% 
  spread(key=Sex,value=total_panicles) %>% 
  dplyr::rename(
    fem_pan ="F",
    mal_pan = "M",
    tempdorm=ztempdorm,
    tempgrow=ztempgrow,
    pptdorm=zpptdorm,
    pptgrow=zpptgrow) %>% 
  subset( fem_pan > 0 & mal_pan > 0 )%>%
  mutate(
    tot_pan=fem_pan+mal_pan,
    osr = fem_pan/tot_pan,
    year_index = year-2014)

#OSR pooled across years
garden_osr_poolyr <- poar.clim_seasonal %>% 
  dplyr::select(site, Sex, zpptgrow,zpptdorm,ztempgrow, ztempdorm,flowerN_t1) %>% 
  group_by(zpptgrow,zpptdorm,ztempgrow,ztempdorm, site, Sex) %>% 
  summarise(total_panicles = sum(flowerN_t1,na.rm=T)) %>% 
  spread(key=Sex,value=total_panicles) %>% 
  rename(fem_pan = "F",
         mal_pan = "M",
         tempdorm=ztempdorm,
         tempgrow=ztempgrow,
         pptdorm=zpptdorm,
         pptgrow=zpptgrow) %>% 
  subset( fem_pan > 0 & mal_pan > 0 )%>%
  mutate(tot_pan=fem_pan+mal_pan,
         osr = fem_pan/tot_pan)

##sex ratio by year
garden_sr <- poar.clim_seasonal %>% 
  dplyr::select( year, site, Sex, zpptgrow,zpptdorm,ztempgrow, ztempdorm,surv_t1) %>% 
  group_by(year, zpptgrow,zpptdorm,ztempgrow, ztempdorm, site, Sex) %>% 
  summarise(total_plants = sum(surv_t1,na.rm=T)) %>% 
  spread(key=Sex,value=total_plants) %>% 
  rename(fem = `F`,
         mal = `M`,
         tempdorm=ztempdorm,
         tempgrow=ztempgrow,
         pptdorm=zpptdorm,
         pptgrow=zpptgrow) %>% 
  subset( fem > 0 & mal > 0 )%>%
  mutate(total=fem+mal,
         sr = fem/total)
##sex ratio pooled across years
garden_sr_poolyr <- poar.clim_seasonal %>%
  dplyr::select(site, Sex, zpptgrow, zpptdorm, ztempgrow, ztempdorm, surv_t1) %>% 
  group_by(zpptgrow,zpptdorm,ztempgrow, ztempdorm, site, Sex) %>% 
  summarise(total_plants = sum(surv_t1,na.rm=T)) %>% 
  spread(key=Sex,value=total_plants) %>% 
  rename(fem = `F`,
         mal = `M`,
         tempdorm=ztempdorm,
         tempgrow=ztempgrow,
         pptdorm=zpptdorm,
         pptgrow=zpptgrow) %>% 
  subset( fem > 0 & mal > 0 )%>%
  mutate(total=fem+mal,
         sr = fem/total)

#mean(garden_osr_poolyr$osr)
garden_osr_dat <- list(y = garden_osr_poolyr$fem_pan,
                       n_trials = garden_osr_poolyr$tot_pan,
                       tempgrow = as.vector(garden_osr_poolyr$tempgrow),
                       tempdorm = as.vector(garden_osr_poolyr$tempdorm),
                       pptdorm=as.vector(garden_osr_poolyr$pptgrow),
                       pptgrow=as.vector(garden_osr_poolyr$pptdorm),
                       N = nrow(garden_osr_poolyr))
# sim_pars <- list(
#   warmup = 1000, 
#   iter = 9000, 
#   thin = 3, 
#   chains = 3
# )

# library(cmdstanr)
# fit_garden_osr <- stan(
#   file = '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Analysis/stan/poar_garden_climate.stan',
#   data = garden_osr_dat,
#   warmup = sim_pars$warmup,
#   iter = sim_pars$iter,
#   thin = sim_pars$thin,
#   chains = sim_pars$chains )

# saveRDS(fit_garden_osr, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Forecasting Models output/fit_garden_osr.rds')
fit_garden_osr <- readRDS(url("https://www.dropbox.com/scl/fi/leogofl0m0gp35xv8zeh1/fit_garden_osr.rds?rlkey=i0b3exnhscocbh14bw7u0kxwe&dl=1"))

traceplot(fit_garden_osr, pars = quote_bare(b0,b_tempgrow,b_pptgrow,b_pptdorm,b_tempgrow,b_tempdorm,b_pptgrow2,b_pptdorm2,b_tempgrow2,b_tempdorm2,b_tempdormpptdorm,b_tempgrowpptgrow))+theme_bw()

## pull out coefficients
coef_garden_osr <- rstan::extract(fit_garden_osr,pars = quote_bare(b0,b_tempgrow,b_pptgrow,b_pptdorm,b_tempdorm,b_pptgrow2,b_pptdorm2,b_tempgrow2,b_tempdorm2,b_tempdormpptdorm,b_tempgrowpptgrow))
#names(coef_garden_osr)

summary(fit_garden_osr,pars = quote_bare(b0,b_tempgrow,b_pptgrow,b_pptdorm,b_tempgrow,b_tempdorm,b_pptgrow2,b_pptdorm2,b_tempgrow2,b_tempdorm2,b_tempdormpptdorm,b_tempgrowpptgrow),probs=c(0.025,0.975))

posterior_osr <- as.array(fit_garden_osr)
color_scheme_set("orange")
osr<-mcmc_intervals(posterior_osr, pars = quote_bare(b0,b_tempgrow,b_pptgrow,b_pptdorm,b_tempdorm,b_pptgrow2,b_pptdorm2,b_tempgrow2,b_tempdorm2,b_tempdormpptdorm,b_tempgrowpptgrow)) + 
  ggplot2::scale_y_discrete(limits = c("b0","b_tempgrow","b_pptgrow","b_pptdorm","b_tempdorm","b_pptgrow2","b_pptdorm2","b_tempgrow2","b_tempdorm2","b_tempdormpptdorm","b_tempgrowpptgrow"),
                            labels=c("b0"="b0",
                                     "b_pptgrow"="pptgrow",
                                     "b_tempgrow"="tempgrow",
                                     "b_pptdorm"="pptdorm",
                                     "b_tempdorm"="tempdorm",
                                     "b_tempdormpptdorm"="tempdormpptdorm",
                                     "b_tempgrowpptgrow"="tempgrowpptgrow",
                                     "b_pptgrow2"=expression(pptgrow^2),
                                     "b_tempgrow2"=expression(tempgrow^2),
                                     "b_pptdorm2"=expression(pptdorm^2),
                                     "b_tempdorm2"=expression(tempdorm^2)))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates ")+
  xlim(-2,2)+
  # ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Posterior_ORS.pdf",useDingbats = F,height=4,width=5)
ggarrange(osr + rremove("ylab"), ncol = 1, nrow = 1)
dev.off()

set.seed(13)
n_post_draws<-300
osr_post_draws <- sample.int(length(coef_garden_osr$b0), n_post_draws)

x_tempgrow <- seq(min(garden_osr_poolyr$tempgrow),max(garden_osr_poolyr$tempgrow),length = 100)
x_pptgrow <- seq(min(garden_osr_poolyr$pptgrow),max(garden_osr_poolyr$pptgrow),length = 100)
x_tempdorm <- seq(min(garden_osr_poolyr$tempdorm),max(garden_osr_poolyr$tempdorm),length = 100)
x_pptdorm <- seq(min(garden_osr_poolyr$pptdorm),max(garden_osr_poolyr$pptdorm),length = 100)

garden_osr_poolyr_plot <- poar.clim_seasonal %>% 
  dplyr::select(site, Sex, pptgrow,pptdorm,tempgrow, tempdorm,flowerN_t1) %>% 
  group_by(pptgrow,pptdorm,tempgrow,tempdorm, site, Sex) %>% 
  summarise(total_panicles = sum(flowerN_t1,na.rm=T)) %>% 
  spread(key=Sex,value=total_panicles) %>% 
  rename(fem_pan = "F",
         mal_pan = "M",
         tempdorm=tempdorm,
         tempgrow=tempgrow,
         pptdorm=pptdorm,
         pptgrow=pptgrow) %>% 
  subset( fem_pan > 0 & mal_pan > 0 )%>%
  mutate(tot_pan=fem_pan+mal_pan,
         osr = fem_pan/tot_pan)

view(garden_osr_poolyr_plot)

mean_coef_osr<- lapply(rstan::extract(fit_garden_osr, pars = quote_bare(b0,b_tempgrow,b_pptgrow,b_pptdorm,b_tempdorm,b_pptgrow2,b_pptdorm2,b_tempgrow2,b_tempdorm2,b_tempdormpptdorm,b_tempgrowpptgrow)) ,mean)



pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/gardens_OSR.pdf",useDingbats = F,width=8,height=4)
par(mar=c(5,5,1.5,0.5),mfrow=c(1,2))

# plot(garden_osr_poolyr$pptgrow*sd(garden_osr_poolyr_plot$pptgrow)+mean(garden_osr_poolyr_plot$pptgrow),garden_osr_poolyr$osr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
#      xlab="Precipitation (grow.season)",ylab="Operational Sex Ratio")
# points(garden_osr_poolyr$pptgrow*sd(garden_osr_poolyr_plot$pptgrow)+mean(garden_osr_poolyr_plot$pptgrow),cex=log(garden_osr_poolyr_plot$tot_pan)*0.5,
#        garden_osr_poolyr$osr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)
# abline(h=0.5,col="gray",lty=2)
# mtext("A",side = 3, adj = 0,cex=1.25)
# for(p in 1:n_post_draws){
#   lines(x_pptgrow*sd(garden_osr_poolyr_plot$pptgrow)+mean(garden_osr_poolyr_plot$pptgrow),
#         invlogit(coef_garden_osr$b0[osr_post_draws[p]] +
#                    coef_garden_osr$b_pptgrow[osr_post_draws[p]] * x_pptgrow +
#                    coef_garden_osr$b_tempgrow[osr_post_draws[p]] * mean(garden_osr_poolyr$tempgrow) +
#                    coef_garden_osr$b_pptdorm[osr_post_draws[p]] * mean(garden_osr_poolyr$pptdorm) + 
#                    coef_garden_osr$b_tempdorm[osr_post_draws[p]] * mean(garden_osr_poolyr$tempdorm) + 
#                    coef_garden_osr$b_pptgrow2[osr_post_draws[p]] * x_pptgrow* x_pptgrow)+
#           coef_garden_osr$b_tempgrow2[osr_post_draws[p]] * mean(garden_osr_poolyr$tempgrow)* mean(garden_osr_poolyr$tempgrow)+
#           coef_garden_osr$b_pptdorm2[osr_post_draws[p]] * mean(garden_osr_poolyr$pptdorm)* mean(garden_osr_poolyr$pptdorm)+
#           coef_garden_osr$b_tempdorm2[osr_post_draws[p]] * mean(garden_osr_poolyr$tempdorm)* mean(garden_osr_poolyr$tempdorm),
#         col=alpha("#E6AB02",0.1))
# }

plot(garden_osr_poolyr$tempgrow*sd(garden_osr_poolyr_plot$tempgrow)+mean(garden_osr_poolyr_plot$tempgrow),garden_osr_poolyr$osr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
     xlab="Temperature (grow.season)",ylab="Operational Sex Ratio")
points(garden_osr_poolyr$tempgrow*sd(garden_osr_poolyr_plot$tempgrow)+mean(garden_osr_poolyr_plot$tempgrow),cex=log(garden_osr_poolyr_plot$tot_pan)*0.5,
       garden_osr_poolyr$osr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)
abline(h=0.5,col="gray",lty=2)
mtext("A",side = 3, adj = 0,cex=1.25)
for(p in 1:n_post_draws){
  lines(x_tempgrow*sd(garden_osr_poolyr_plot$tempgrow)+mean(garden_osr_poolyr_plot$tempgrow),
        invlogit(coef_garden_osr$b0[osr_post_draws[p]] +
                   coef_garden_osr$b_tempgrow[osr_post_draws[p]] * x_tempgrow +
                   coef_garden_osr$b_pptgrow[osr_post_draws[p]] * mean(garden_osr_poolyr$pptgrow) +
                   coef_garden_osr$b_pptdorm[osr_post_draws[p]] * mean(garden_osr_poolyr$pptdorm) + 
                   coef_garden_osr$b_tempdorm[osr_post_draws[p]] * mean(garden_osr_poolyr$tempdorm) + 
                   coef_garden_osr$b_pptgrow2[osr_post_draws[p]] * mean(garden_osr_poolyr$pptgrow)* mean(garden_osr_poolyr$pptgrow))+
          coef_garden_osr$b_tempgrow2[osr_post_draws[p]] * x_tempgrow* x_tempgrow+
          coef_garden_osr$b_pptdorm2[osr_post_draws[p]] * mean(garden_osr_poolyr$pptdorm)* mean(garden_osr_poolyr$pptdorm)+
          coef_garden_osr$b_tempdorm2[osr_post_draws[p]] * mean(garden_osr_poolyr$tempdorm)* mean(garden_osr_poolyr$tempdorm),
        col=alpha("#E6AB02",0.1))
}


# plot(garden_osr_poolyr$pptdorm*sd(garden_osr_poolyr_plot$pptdorm)+mean(garden_osr_poolyr_plot$pptdorm),garden_osr_poolyr$osr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
#      xlab="Precipitation (dorm.season)",ylab="Operational Sex Ratio")
# points(garden_osr_poolyr$pptdorm*sd(garden_osr_poolyr_plot$pptdorm)+mean(garden_osr_poolyr_plot$pptdorm),cex=log(garden_osr_poolyr_plot$tot_pan)*0.5,
#        garden_osr_poolyr$osr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)
# 
# mtext("C",side = 3, adj = 0,cex=1.25)
# abline(h=0.5,col="gray",lty=2)
# for(p in 1:n_post_draws){
#   lines(x_pptdorm*sd(garden_osr_poolyr_plot$pptdorm)+mean(garden_osr_poolyr_plot$pptdorm),
#         invlogit(coef_garden_osr$b0[osr_post_draws[p]] +
#                    coef_garden_osr$b_pptgrow[osr_post_draws[p]] * mean(garden_osr_poolyr$pptgrow) +
#                    coef_garden_osr$b_tempgrow[osr_post_draws[p]] * mean(garden_osr_poolyr$tempgrow) +
#                    coef_garden_osr$b_pptdorm[osr_post_draws[p]] * x_pptdorm + 
#                    coef_garden_osr$b_tempdorm[osr_post_draws[p]] * mean(garden_osr_poolyr$tempdorm) + 
#                    coef_garden_osr$b_pptgrow2[osr_post_draws[p]] * mean(garden_osr_poolyr$pptgrow)* mean(garden_osr_poolyr$pptgrow))+
#           coef_garden_osr$b_tempgrow2[osr_post_draws[p]] * mean(garden_osr_poolyr$tempgrow)* mean(garden_osr_poolyr$tempgrow)+
#           coef_garden_osr$b_pptdorm2[osr_post_draws[p]] * x_pptdorm* x_pptdorm+
#           coef_garden_osr$b_tempdorm2[osr_post_draws[p]] * mean(garden_osr_poolyr$tempdorm)* mean(garden_osr_poolyr$tempdorm),
#         col=alpha("#E6AB02",0.1))
# }

plot(garden_osr_poolyr$tempdorm*sd(garden_osr_poolyr_plot$tempdorm)+mean(garden_osr_poolyr_plot$tempdorm),garden_osr_poolyr$osr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
     xlab="Temperature (dorm.season)",ylab="Operational Sex Ratio")
points(garden_osr_poolyr$tempdorm*sd(garden_osr_poolyr_plot$tempdorm)+mean(garden_osr_poolyr_plot$tempdorm),cex=log(garden_osr_poolyr_plot$tot_pan)*0.5,
       garden_osr_poolyr$osr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)
mtext("B",side = 3, adj = 0,cex=1.25)
abline(h=0.5,col="gray",lty=2)
for(p in 1:n_post_draws){
  lines(x_tempdorm*sd(garden_osr_poolyr_plot$tempdorm)+mean(garden_osr_poolyr_plot$tempdorm),
        invlogit(coef_garden_osr$b0[osr_post_draws[p]] +
                   coef_garden_osr$b_tempdorm[osr_post_draws[p]] * x_tempdorm +
                   coef_garden_osr$b_tempgrow[osr_post_draws[p]] * mean(garden_osr_poolyr$tempgrow) +
                   coef_garden_osr$b_pptdorm[osr_post_draws[p]] * mean(garden_osr_poolyr$pptdorm) + 
                   coef_garden_osr$b_pptgrow[osr_post_draws[p]] * mean(garden_osr_poolyr$pptgrow) + 
                   coef_garden_osr$b_pptgrow2[osr_post_draws[p]] * mean(garden_osr_poolyr$pptgrow)* mean(garden_osr_poolyr$pptgrow))+
          coef_garden_osr$b_tempgrow2[osr_post_draws[p]] * mean(garden_osr_poolyr$tempgrow)* mean(garden_osr_poolyr$tempgrow)+
          coef_garden_osr$b_pptdorm2[osr_post_draws[p]] * mean(garden_osr_poolyr$pptdorm)* mean(garden_osr_poolyr$pptdorm)+
          coef_garden_osr$b_tempdorm2[osr_post_draws[p]] * x_tempdorm* x_tempdorm,
        col=alpha("#E6AB02",0.1))
}


dev.off()


garden_sr_dat <- list(y = garden_sr_poolyr$fem,
                      n_trials = garden_sr_poolyr$total,
                      tempgrow = as.vector(garden_osr_poolyr$tempgrow),
                      tempdorm = as.vector(garden_osr_poolyr$tempdorm),
                      pptgrow = as.vector(garden_osr_poolyr$pptgrow),
                      pptdorm = as.vector(garden_osr_poolyr$pptdorm),
                      N = nrow(garden_sr_poolyr))


pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/gardens_OSR_ESA.pdf",useDingbats = F,width=8,height=4)
par(mar=c(5,5,1.5,0.5),mfrow=c(1,2))

# plot(garden_osr_poolyr$pptgrow*sd(garden_osr_poolyr_plot$pptgrow)+mean(garden_osr_poolyr_plot$pptgrow),garden_osr_poolyr$osr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
#      xlab="Precipitation (grow.season)",ylab="Sex ratio")
# abline(h=0.5,col="gray",lty=2)
# mtext("A",side = 3, adj = 0,cex=1.25)
# for(p in 1:n_post_draws){
#   lines(x_pptgrow*sd(garden_osr_poolyr_plot$pptgrow)+mean(garden_osr_poolyr_plot$pptgrow),
#         invlogit(coef_garden_osr_pptgrow$b0[osr_post_draws[p]] +
#                    coef_garden_osr_pptgrow$b_pptgrow[osr_post_draws[p]] * x_pptgrow +
#                    coef_garden_osr_pptgrow$b_pptgrow2[osr_post_draws[p]] * x_pptgrow* x_pptgrow),
#         col=alpha("#E6AB02",0.1))
# }
# points(garden_osr_poolyr$pptgrow*sd(garden_osr_poolyr_plot$pptgrow)+mean(garden_osr_poolyr_plot$pptgrow),cex=1.5,
#        garden_osr_poolyr$osr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)
# 


plot(garden_osr_poolyr$tempgrow*sd(garden_osr_poolyr_plot$tempgrow)+mean(garden_osr_poolyr_plot$tempgrow),garden_osr_poolyr$osr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
     xlab="Temperature (grow.season)",ylab="Sex ratio (proportion of female)")
abline(h=0.5,col="gray",lty=2)
mtext("A",side = 3, adj = 0,cex=1.25)
for(p in 1:n_post_draws){
  lines(x_tempgrow*sd(garden_osr_poolyr_plot$tempgrow)+mean(garden_osr_poolyr_plot$tempgrow),
        invlogit(coef_garden_osr_tempgrow$b0[osr_post_draws[p]] +
                   coef_garden_osr_tempgrow$b_tempgrow[osr_post_draws[p]] * x_tempgrow +
                   coef_garden_osr_tempgrow$b_tempgrow2[osr_post_draws[p]] * x_tempgrow* x_tempgrow),
        col=alpha("#E6AB02",0.1))
}
points(garden_osr_poolyr$tempgrow*sd(garden_osr_poolyr_plot$tempgrow)+mean(garden_osr_poolyr_plot$tempgrow),cex=1.5,
       garden_osr_poolyr$osr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)



# plot(garden_osr_poolyr$pptdorm*sd(garden_osr_poolyr_plot$pptdorm)+mean(garden_osr_poolyr_plot$pptdorm),garden_osr_poolyr$osr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
#      xlab="Precipitation (dorm.season)",ylab="Sex ratio")
# mtext("C",side = 3, adj = 0,cex=1.25)
# abline(h=0.5,col="gray",lty=2)
# for(p in 1:n_post_draws){
#   lines(x_pptdorm*sd(garden_osr_poolyr_plot$pptdorm)+mean(garden_osr_poolyr_plot$pptdorm),
#         invlogit(coef_garden_osr_pptdorm$b0[osr_post_draws[p]] +
#                    coef_garden_osr_pptdorm$b_pptdorm[osr_post_draws[p]] * x_pptdorm +
#                    coef_garden_osr_pptdorm$b_pptdorm2[osr_post_draws[p]] * x_pptdorm* x_pptdorm),
#         col=alpha("#E6AB02",0.1))
# }
# points(garden_osr_poolyr$pptdorm*sd(garden_osr_poolyr_plot$pptdorm)+mean(garden_osr_poolyr_plot$pptdorm),cex=1.5,
#        garden_osr_poolyr$osr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)


plot(garden_osr_poolyr$tempdorm*sd(garden_osr_poolyr_plot$tempdorm)+mean(garden_osr_poolyr_plot$tempdorm),garden_osr_poolyr$osr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
     xlab="Temperature (dorm.season)",ylab="Sex ratio (proportion of female)")
mtext("B",side = 3, adj = 0,cex=1.25)
abline(h=0.5,col="gray",lty=2)
for(p in 1:n_post_draws){
  lines(x_tempdorm*sd(garden_osr_poolyr_plot$tempdorm)+mean(garden_osr_poolyr_plot$tempdorm),
        invlogit(coef_garden_osr_tempdorm$b0[osr_post_draws[p]] +
                   coef_garden_osr_tempdorm$b_tempdorm[osr_post_draws[p]] * x_tempdorm +
                   coef_garden_osr_tempdorm$b_tempdorm2[osr_post_draws[p]] * x_tempdorm* x_tempdorm),
        col=alpha("#E6AB02",0.1))
}
points(garden_osr_poolyr$tempdorm*sd(garden_osr_poolyr_plot$tempdorm)+mean(garden_osr_poolyr_plot$tempdorm),cex=1.5,
       garden_osr_poolyr$osr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)

dev.off()


garden_sr_dat <- list(y = garden_sr_poolyr$fem,
                      n_trials = garden_sr_poolyr$total,
                      tempgrow = as.vector(garden_sr_poolyr$tempgrow),
                      tempdorm = as.vector(garden_sr_poolyr$tempdorm),
                      pptgrow = as.vector(garden_sr_poolyr$pptgrow),
                      pptdorm = as.vector(garden_sr_poolyr$pptdorm),
                      N = nrow(garden_sr_poolyr))


fit_garden_sr <- stan(
  file = '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Analysis/stan/poar_garden_climate.stan',
  data = garden_sr_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )

saveRDS(fit_garden_sr, '/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Forecasting Models output/fit_garden_sr.rds')
fit_garden_sr <- readRDS(url("https://www.dropbox.com/scl/fi/xplzsxmxxwx2scqjzblom/fit_garden_sr.rds?rlkey=t514ak5c605egen4f1v4bmom9&dl=1"))

traceplot(fit_garden_sr, pars = quote_bare(b0,b_tempgrow,b_pptgrow,b_pptdorm,b_tempgrow,b_tempdorm,b_pptgrow2,b_pptdorm2,b_tempgrow2,b_tempdorm2,b_tempdormpptdorm,b_tempgrowpptgrow))+theme_bw()

## pull out coefficients
coef_garden_sr <- rstan::extract(fit_garden_sr,pars = quote_bare(b0,b_tempgrow,b_pptgrow,b_pptdorm,b_tempdorm,b_pptgrow2,b_pptdorm2,b_tempgrow2,b_tempdorm2,b_tempdormpptdorm,b_tempgrowpptgrow))

posterior_sr <- as.array(fit_garden_sr)
color_scheme_set("orange")
sr<-mcmc_intervals(posterior_sr, pars = quote_bare(b0,b_tempgrow,b_pptgrow,b_pptdorm,b_tempdorm,b_pptgrow2,b_pptdorm2,b_tempgrow2,b_tempdorm2,b_tempdormpptdorm,b_tempgrowpptgrow)) + 
  ggplot2::scale_y_discrete(limits = c("b0","b_tempgrow","b_pptgrow","b_pptdorm","b_tempdorm","b_pptgrow2","b_pptdorm2","b_tempgrow2","b_tempdorm2","b_tempdormpptdorm","b_tempgrowpptgrow"),
                            labels=c("b0"="b0",
                                     "b_pptgrow"="pptgrow",
                                     "b_tempgrow"="tempgrow",
                                     "b_pptdorm"="pptdorm",
                                     "b_tempdorm"="tempdorm",
                                     "b_tempdormpptdorm"="tempdormpptdorm",
                                     "b_tempgrowpptgrow"="tempgrowpptgrow",
                                     "b_pptgrow2"=expression(pptgrow^2),
                                     "b_tempgrow2"=expression(tempgrow^2),
                                     "b_pptdorm2"=expression(pptdorm^2),
                                     "b_tempdorm2"=expression(tempdorm^2)))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates ")+
  xlim(-1,1)+
  # ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))


pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Posterior_SR.pdf",useDingbats = F,height=4,width=5)
ggarrange(sr + rremove("ylab"), ncol = 1, nrow = 1)
dev.off()

set.seed(13)
n_post_draws<-300
sr_post_draws <- sample.int(length(coef_garden_sr$b0), n_post_draws)

garden_sr_poolyr_plot <- poar.clim_seasonal %>%
  dplyr::select(site, Sex, pptgrow, pptdorm, tempgrow, tempdorm, surv_t1) %>% 
  group_by(pptgrow,pptdorm,tempgrow, tempdorm, site, Sex) %>% 
  summarise(total_plants = sum(surv_t1,na.rm=T)) %>% 
  spread(key=Sex,value=total_plants) %>% 
  rename(fem = `F`,
         mal = `M`,
         tempdorm=tempdorm,
         tempgrow=tempgrow,
         pptdorm=pptdorm,
         pptgrow=pptgrow) %>% 
  subset( fem > 0 & mal > 0 )%>%
  mutate(total=fem+mal,
         sr = fem/total)


pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/gardens_SR.pdf",useDingbats = F,width=9,height=8)
par(mar=c(5,5,1.5,0.5),mfrow=c(2,2))

plot(garden_sr_poolyr$pptgrow*sd(garden_sr_poolyr_plot$pptgrow)+mean(garden_sr_poolyr_plot$pptgrow),garden_sr_poolyr$sr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
     xlab="Precipitation (grow.season)",ylab="Proportion of female plants")
points(garden_sr_poolyr$pptgrow*sd(garden_sr_poolyr_plot$pptgrow)+mean(garden_sr_poolyr_plot$pptgrow),cex=log(garden_sr_poolyr_plot$total)*0.5,
       garden_sr_poolyr$sr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)
abline(h=0.5,col="gray",lty=2)
mtext("A",side = 3, adj = 0,cex=1.25)
for(p in 1:n_post_draws){
  lines(x_pptgrow*sd(garden_sr_poolyr_plot$pptgrow)+mean(garden_sr_poolyr_plot$pptgrow),
        invlogit(coef_garden_sr$b0[osr_post_draws[p]] +
                   coef_garden_sr$b_pptgrow[sr_post_draws[p]] * x_pptgrow +
                   coef_garden_sr$b_tempgrow[sr_post_draws[p]] * mean(garden_sr_poolyr$tempgrow) +
                   coef_garden_sr$b_pptdorm[sr_post_draws[p]] * mean(garden_sr_poolyr$pptdorm) + 
                   coef_garden_sr$b_tempdorm[sr_post_draws[p]] * mean(garden_sr_poolyr$tempdorm) + 
                   coef_garden_sr$b_pptgrow2[sr_post_draws[p]] * x_pptgrow* x_pptgrow)+
          coef_garden_osr$b_tempgrow2[sr_post_draws[p]] * mean(garden_sr_poolyr$tempgrow)* mean(garden_sr_poolyr$tempgrow)+
          coef_garden_sr$b_pptdorm2[sr_post_draws[p]] * mean(garden_sr_poolyr$pptdorm)* mean(garden_sr_poolyr$pptdorm)+
          coef_garden_sr$b_tempdorm2[sr_post_draws[p]] * mean(garden_sr_poolyr$tempdorm)* mean(garden_sr_poolyr$tempdorm),
        col=alpha("#E6AB02",0.1))
}



plot(garden_sr_poolyr$tempgrow*sd(garden_sr_poolyr_plot$tempgrow)+mean(garden_sr_poolyr_plot$tempgrow),garden_sr_poolyr$sr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
     xlab="Temerature (grow.season)",ylab="Proportion of female plants")
abline(h=0.5,col="gray",lty=2)
mtext("B",side = 3, adj = 0,cex=1.25)
points(garden_sr_poolyr$tempgrow*sd(garden_sr_poolyr_plot$tempgrow)+mean(garden_sr_poolyr_plot$tempgrow),cex=log(garden_sr_poolyr_plot$total)*0.5,
       garden_sr_poolyr$sr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)
for(p in 1:n_post_draws){
  lines(x_tempgrow*sd(garden_sr_poolyr_plot$tempgrow)+mean(garden_sr_poolyr_plot$tempgrow),
        invlogit(coef_garden_sr$b0[sr_post_draws[p]] +
                   coef_garden_sr$b_tempgrow[sr_post_draws[p]] * x_tempgrow +
                   coef_garden_sr$b_pptgrow[sr_post_draws[p]] * mean(garden_sr_poolyr$pptgrow) +
                   coef_garden_sr$b_pptdorm[sr_post_draws[p]] * mean(garden_sr_poolyr$pptdorm) + 
                   coef_garden_sr$b_tempdorm[sr_post_draws[p]] * mean(garden_sr_poolyr$tempdorm) + 
                   coef_garden_sr$b_pptgrow2[sr_post_draws[p]] * mean(garden_sr_poolyr$pptgrow)* mean(garden_sr_poolyr$pptgrow))+
          coef_garden_sr$b_tempgrow2[sr_post_draws[p]] * x_tempgrow* x_tempgrow+
          coef_garden_sr$b_pptdorm2[sr_post_draws[p]] * mean(garden_sr_poolyr$pptdorm)* mean(garden_sr_poolyr$pptdorm)+
          coef_garden_sr$b_tempdorm2[sr_post_draws[p]] * mean(garden_sr_poolyr$tempdorm)* mean(garden_sr_poolyr$tempdorm),
        col=alpha("#E6AB02",0.1))
}



plot(garden_sr_poolyr$pptdorm*sd(garden_sr_poolyr_plot$pptdorm)+mean(garden_sr_poolyr_plot$pptdorm),garden_sr_poolyr$sr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
     xlab="Precipitation (dorm.season)",ylab="Proportion of female plants")
mtext("C",side = 3, adj = 0,cex=1.25)
abline(h=0.5,col="gray",lty=2)
points(garden_sr_poolyr$pptdorm*sd(garden_sr_poolyr_plot$pptdorm)+mean(garden_sr_poolyr_plot$pptdorm),cex=log(garden_sr_poolyr_plot$total)*0.5,
       garden_sr_poolyr$sr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)
for(p in 1:n_post_draws){
  lines(x_pptdorm*sd(garden_sr_poolyr_plot$pptdorm)+mean(garden_sr_poolyr_plot$pptdorm),
        invlogit(coef_garden_sr$b0[osr_post_draws[p]] +
                   coef_garden_sr$b_pptdorm[sr_post_draws[p]] * x_pptdorm +
                   coef_garden_sr$b_tempgrow[sr_post_draws[p]] * mean(garden_sr_poolyr$tempgrow) +
                   coef_garden_sr$b_pptdorm[sr_post_draws[p]] * mean(garden_sr_poolyr$pptdorm) + 
                   coef_garden_sr$b_tempdorm[sr_post_draws[p]] * mean(garden_sr_poolyr$tempdorm) + 
                   coef_garden_sr$b_pptgrow2[sr_post_draws[p]] * mean(garden_sr_poolyr$pptgrow)* mean(garden_sr_poolyr$pptgrow))+
          coef_garden_sr$b_tempgrow2[sr_post_draws[p]] * mean(garden_sr_poolyr$tempgrow)* mean(garden_sr_poolyr$tempgrow)+
          coef_garden_sr$b_pptdorm2[sr_post_draws[p]] * mean(garden_sr_poolyr$pptdorm)* mean(garden_sr_poolyr$pptdorm)+
          coef_garden_sr$b_tempdorm2[sr_post_draws[p]] * mean(garden_sr_poolyr$tempdorm)* mean(garden_sr_poolyr$tempdorm),
        col=alpha("#E6AB02",0.1))
}


plot(garden_sr_poolyr$tempdorm*sd(garden_sr_poolyr_plot$tempdorm)+mean(garden_sr_poolyr_plot$tempdorm),garden_sr_poolyr$sr,ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,
     xlab="Temperature (dorm.season)",ylab="Proportion of female plants")
mtext("D",side = 3, adj = 0,cex=1.25)
abline(h=0.5,col="gray",lty=2)
points(garden_sr_poolyr$tempdorm*sd(garden_sr_poolyr_plot$tempdorm)+mean(garden_sr_poolyr_plot$tempdorm),cex=log(garden_sr_poolyr_plot$total)*0.5,
       garden_sr_poolyr$sr,lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)
for(p in 1:n_post_draws){
  lines(x_tempdorm*sd(garden_sr_poolyr_plot$tempdorm)+mean(garden_sr_poolyr_plot$tempdorm),
        invlogit(coef_garden_sr$b0[sr_post_draws[p]] +
                   coef_garden_sr$b_tempdorm[sr_post_draws[p]] * x_tempdorm +
                   coef_garden_sr$b_tempgrow[sr_post_draws[p]] * mean(garden_sr_poolyr$tempgrow) +
                   coef_garden_sr$b_pptdorm[sr_post_draws[p]] * mean(garden_sr_poolyr$pptdorm) + 
                   coef_garden_sr$b_pptgrow[sr_post_draws[p]] * mean(garden_sr_poolyr$pptgrow) + 
                   coef_garden_sr$b_pptgrow2[sr_post_draws[p]] * mean(garden_sr_poolyr$pptgrow)* mean(garden_osr_poolyr$pptgrow))+
          coef_garden_sr$b_tempgrow2[sr_post_draws[p]] * mean(garden_sr_poolyr$tempgrow)* mean(garden_sr_poolyr$tempgrow)+
          coef_garden_sr$b_pptdorm2[sr_post_draws[p]] * mean(garden_sr_poolyr$pptdorm)* mean(garden_sr_poolyr$pptdorm)+
          coef_garden_sr$b_tempdorm2[sr_post_draws[p]] * x_tempdorm* x_tempdorm,
        col=alpha("#E6AB02",0.1))
}

dev.off()


# Vital rates informatiom ----
fit_allsites_season <- readRDS(url("https://www.dropbox.com/scl/fi/1sbhjorfv15ianvdwli96/fit_1.5_0.5.rds?rlkey=rm9sxy29qzc7himps7n6y0xw6&dl=1")) #ignore warning from readRDS

# pull out stan coefficients
poar <- poar.clim_seasonal  
fit_full <- fit_allsites_season  
surv_coef <- rstan::extract(fit_full, pars = quote_bare(b0_s,bsizesex_s,btempdormpptdorm_s,btempgrowpptgrow_s,bsex_s,
                                                        bsize_s,bpptgrow_s,bpptdorm_s,btempgrow_s,btempdorm_s,
                                                        bpptgrowsex_s,bpptdormsex_s,btempgrowsex_s,btempdormsex_s,
                                                        btempdormpptdormsex_s,btempgrowpptgrowsex_s,bpptgrow2_s,
                                                        bpptdorm2_s,btempgrow2_s,btempdorm2_s,bpptgrow2sex_s,bpptdorm2sex_s,
                                                        btempgrow2sex_s,btempdorm2sex_s,
                                                        site_tau_s,block_tau_s,source_tau_s))

grow_coef <- rstan::extract(fit_full, pars = quote_bare(b0_g,btempdormpptdorm_g,btempgrowpptgrow_g,bsizesex_g,bsex_g,
                                                        bsize_g,bpptgrow_g,bpptdorm_g,btempgrow_g,btempdorm_g,bpptgrowsex_g,
                                                        bpptdormsex_g,btempgrowsex_g,btempdormsex_g,btempdormpptdormsex_g,
                                                        btempgrowpptgrowsex_g,bpptgrow2_g,bpptdorm2_g,btempgrow2_g,btempdorm2_g,
                                                        bpptgrow2sex_g,bpptdorm2sex_g,btempgrow2sex_g,btempdorm2sex_g,
                                                        site_tau_g,block_tau_g,source_tau_g,sigma))

flow_coef <- rstan::extract(fit_full, pars = quote_bare(b0_f,bsizesex_f,bsex_f, bsize_f,bpptgrow_f,bpptdorm_f,btempgrow_f,
                                                        btempdorm_f,bpptgrowsex_f,bpptdormsex_f,btempgrowsex_f,btempdormsex_f,
                                                        btempdormpptdorm_f,btempgrowpptgrow_f,btempdormpptdormsex_f,
                                                        btempgrowpptgrowsex_f,bpptgrow2_f,bpptdorm2_f,btempgrow2_f,
                                                        btempdorm2_f,bpptgrow2sex_f,bpptdorm2sex_f,btempgrow2sex_f,btempdorm2sex_f,
                                                        site_tau_f,block_tau_f,source_tau_f))

panic_coef <- rstan::extract(fit_full, pars = quote_bare(b0_p,bsex_p,bsex_p,bsizesex_p, bsize_p,bpptgrow_p,bpptdorm_p,
                                                         btempgrow_p,btempdorm_p,bpptgrowsex_p,bpptdormsex_p,btempgrowsex_p,
                                                         btempdormsex_p,btempdormpptdorm_p,btempgrowpptgrow_p,btempdormpptdormsex_p,
                                                         btempgrowpptgrowsex_p,bpptgrow2_p,bpptdorm2_p,btempgrow2_p,
                                                         btempdorm2_p,bpptgrow2sex_p,bpptdorm2sex_p,btempgrow2sex_p,btempdorm2sex_p,
                                                         site_tau_p,block_tau_p,source_tau_p))

via_coef <- rstan::extract(fit_full, pars = quote_bare(v0,a_v,m,lambda_d))

# Matrix model (lambda vs climate)----
##load MPM functions
source("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Analysis/twosexMPMLTRE.R")
# source("/Users/tm9/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Analysis/twosexMPMLTRE.R")

## read in the seedling survival data from Poa autumnalis
POAU <- read.csv("https://www.dropbox.com/s/7f2yfz49fzrcp5i/POAU.csv?dl=1")
sdlg_surv <- POAU %>% filter(year_recruit==year_t) %>% summarise(sdlg_surv = mean(spring_survival_t1,na.rm=T))

## pull out mean stan coefficients
mean_coefs<- lapply(rstan::extract(fit_full, pars = quote_bare(b0_s,bsizesex_s,btempdormpptdorm_s,btempgrowpptgrow_s,bsex_s,
                                                               bsize_s,bpptgrow_s,bpptdorm_s,btempgrow_s,btempdorm_s,
                                                               bpptgrowsex_s,bpptdormsex_s,btempgrowsex_s,btempdormsex_s,
                                                               btempdormpptdormsex_s,btempgrowpptgrowsex_s,bpptgrow2_s,
                                                               bpptdorm2_s,btempgrow2_s,btempdorm2_s,bpptgrow2sex_s,bpptdorm2sex_s,
                                                               btempgrow2sex_s,btempdorm2sex_s,
                                                               site_tau_s,block_tau_s,source_tau_s)) ,mean)

mean_coefg<- lapply(rstan::extract(fit_full, pars = quote_bare(b0_g,btempdormpptdorm_g,btempgrowpptgrow_g,bsizesex_g,bsex_g,
                                                               bsize_g,bpptgrow_g,bpptdorm_g,btempgrow_g,btempdorm_g,bpptgrowsex_g,
                                                               bpptdormsex_g,btempgrowsex_g,btempdormsex_g,btempdormpptdormsex_g,
                                                               btempgrowpptgrowsex_g,bpptgrow2_g,bpptdorm2_g,btempgrow2_g,
                                                               btempdorm2_g,bpptgrow2sex_g,bpptdorm2sex_g,btempgrow2sex_g,btempdorm2sex_g,
                                                               site_tau_g,block_tau_g,source_tau_g,sigma)) ,mean)

mean_coeff<- lapply(rstan::extract(fit_full, pars = quote_bare(b0_f,bsizesex_f,bsex_f, bsize_f,bpptgrow_f,bpptdorm_f,btempgrow_f,
                                                               btempdorm_f,bpptgrowsex_f,bpptdormsex_f,btempgrowsex_f,btempdormsex_f,
                                                               btempdormpptdorm_f,btempgrowpptgrow_f,btempdormpptdormsex_f,
                                                               btempgrowpptgrowsex_f,bpptgrow2_f,bpptdorm2_f,btempgrow2_f,
                                                               btempdorm2_f,bpptgrow2sex_f,bpptdorm2sex_f,btempgrow2sex_f,btempdorm2sex_f,
                                                               site_tau_f,block_tau_f,source_tau_f)) ,mean)

mean_coefp<- lapply(rstan::extract(fit_full, pars = quote_bare(b0_p,bsex_p,bsex_p,bsizesex_p, bsize_p,bpptgrow_p,bpptdorm_p,
                                                               btempgrow_p,btempdorm_p,bpptgrowsex_p,bpptdormsex_p,btempgrowsex_p,
                                                               btempdormsex_p,btempdormpptdorm_p,btempgrowpptgrow_p,btempdormpptdormsex_p,
                                                               btempgrowpptgrowsex_p,bpptgrow2_p,bpptdorm2_p,btempgrow2_p,
                                                               btempdorm2_p,bpptgrow2sex_p,bpptdorm2sex_p,btempgrow2sex_p,btempdorm2sex_p,
                                                               site_tau_p,block_tau_p,source_tau_p)) ,mean)

mean_coefv<-lapply(rstan::extract(fit_full, pars = quote_bare(v0,a_v,m,lambda_d)) ,mean)

# set up female and male parameter vectors to run the Matrix
F_params <- M_params <- c()
## survival
F_params$surv_mu <- mean_coefs$b0_s
F_params$surv_size <- mean_coefs$bsize_s
F_params$surv_pptgrow <- mean_coefs$bpptgrow_s 
F_params$surv_pptdorm <- mean_coefs$bpptdorm_s 
F_params$surv_tempgrow <- mean_coefs$btempgrow_s 
F_params$surv_tempdorm <- mean_coefs$btempdorm_s  
F_params$surv_tempdorm_pptdorm<-mean_coefs$btempdormpptdorm_s
F_params$surv_tempgrow_pptgrow<-mean_coefs$btempgrowpptgrow_s
F_params$surv_pptgrow2<-mean_coefs$bpptgrow2_s
F_params$surv_pptdorm2<-mean_coefs$bpptdorm2_s
F_params$surv_tempgrow2<-mean_coefs$btempgrow2_s
F_params$surv_tempdorm2<-mean_coefs$btempdorm2_s
M_params$surv_mu <- mean_coefs$b0_s + mean_coefs$bsex_s  
M_params$surv_size <- mean_coefs$bsize_s + mean_coefs$bsizesex_s
M_params$surv_pptgrow <- mean_coefs$bpptgrow_s  + mean_coefs$bpptgrowsex_s 
M_params$surv_pptdorm <- mean_coefs$bpptdorm_s  + mean_coefs$bpptdormsex_s 
M_params$surv_tempgrow <- mean_coefs$btempgrow_s + mean_coefs$btempgrowsex_s #
M_params$surv_tempdorm <-  mean_coefs$btempdorm_s  + mean_coefs$btempdormsex_s  
M_params$surv_tempdorm_pptdorm<-mean_coefs$btempdormpptdorm_s + mean_coefs$btempdormpptdormsex_s
M_params$surv_tempgrow_pptgrow<-mean_coefs$btempgrowpptgrow_s + mean_coefs$btempgrowpptgrowsex_s
M_params$surv_pptgrow2<-mean_coefs$bpptgrow2_s + mean_coefs$bpptgrow2sex_s
M_params$surv_pptdorm2<-mean_coefs$bpptdorm2_s + mean_coefs$bpptdorm2sex_s
M_params$surv_tempgrow2<-mean_coefs$btempgrow2_s + mean_coefs$btempgrow2sex_s
M_params$surv_tempdorm2<-mean_coefs$btempdorm2_s + mean_coefs$btempdorm2sex_s
## growth
F_params$grow_mu <- mean_coefg$b0_g
F_params$grow_size <- mean_coefg$bsize_g
F_params$grow_pptgrow <- mean_coefg$bpptgrow_g 
F_params$grow_pptdorm <- mean_coefg$bpptdorm_g 
F_params$grow_tempgrow <- mean_coefg$btempgrow_g 
F_params$grow_tempdorm <- mean_coefg$btempdorm_g  
F_params$grow_tempdorm_pptdorm<-mean_coefg$btempdormpptdorm_g
F_params$grow_tempgrow_pptgrow<-mean_coefg$btempgrowpptgrow_g
F_params$grow_pptgrow2<-mean_coefg$bpptgrow2_g
F_params$grow_pptdorm2<-mean_coefg$bpptdorm2_g
F_params$grow_tempgrow2<-mean_coefg$btempgrow2_g
F_params$grow_tempdorm2<-mean_coefg$btempdorm2_g
F_params$sigma_g <- mean_coefg$sigma 
M_params$grow_mu <- mean_coefg$b0_g + mean_coefg$bsex_g  
M_params$grow_size <- mean_coefg$bsize_g + mean_coefg$bsizesex_g
M_params$grow_pptgrow <- mean_coefg$bpptgrow_g  + mean_coefg$bpptgrowsex_g 
M_params$grow_pptdorm <- mean_coefg$bpptdorm_g  + mean_coefg$bpptdormsex_g 
M_params$grow_tempgrow <- mean_coefg$btempgrow_g + mean_coefg$btempgrowsex_g #
M_params$grow_tempdorm <-  mean_coefg$btempdorm_g  + mean_coefg$btempdormsex_g  
M_params$grow_tempdorm_pptdorm<-mean_coefg$btempdormpptdorm_g + mean_coefg$btempdormpptdormsex_g
M_params$grow_tempgrow_pptgrow<-mean_coefg$btempgrowpptgrow_g + mean_coefg$btempgrowpptgrowsex_g
M_params$grow_pptgrow2<-mean_coefg$bpptgrow2_g + mean_coefg$bpptgrow2sex_g
M_params$grow_pptdorm2<-mean_coefg$bpptdorm2_g + mean_coefg$bpptdorm2sex_g
M_params$grow_tempgrow2<-mean_coefg$btempgrow2_g + mean_coefg$btempgrow2sex_g
M_params$grow_tempdorm2<-mean_coefg$btempdorm2_g + mean_coefg$btempdorm2sex_g
M_params$sigma_g <- mean_coefg$sigma 
## flowering
F_params$flow_mu <- mean_coeff$b0_f
F_params$flow_size <- mean_coeff$bsize_f
F_params$flow_pptgrow <- mean_coeff$bpptgrow_f 
F_params$flow_pptdorm <- mean_coeff$bpptdorm_f 
F_params$flow_tempgrow <- mean_coeff$btempgrow_f 
F_params$flow_tempdorm <- mean_coeff$btempdorm_f  
F_params$flow_tempdorm_pptdorm<-mean_coeff$btempdormpptdorm_f
F_params$flow_tempgrow_pptgrow<-mean_coeff$btempgrowpptgrow_f
F_params$flow_pptgrow2<-mean_coeff$bpptgrow2_f
F_params$flow_pptdorm2<-mean_coeff$bpptdorm2_f
F_params$flow_tempgrow2<-mean_coeff$btempgrow2_f
F_params$flow_tempdorm2<-mean_coeff$btempdorm2_f
M_params$flow_mu <- mean_coeff$b0_f + mean_coeff$bsex_f  
M_params$flow_size <- mean_coeff$bsize_f + mean_coeff$bsizesex_f
M_params$flow_pptgrow <- mean_coeff$bpptgrow_f  + mean_coeff$bpptgrowsex_f 
M_params$flow_pptdorm <- mean_coeff$bpptdorm_f  + mean_coeff$bpptdormsex_f 
M_params$flow_tempgrow <- mean_coeff$btempgrow_f + mean_coeff$btempgrowsex_f #
M_params$flow_tempdorm <-  mean_coeff$btempdorm_f  + mean_coeff$btempdormsex_f  
M_params$flow_tempdorm_pptdorm<-mean_coeff$btempdormpptdorm_f + mean_coeff$btempdormpptdormsex_f
M_params$flow_tempgrow_pptgrow<-mean_coeff$btempgrowpptgrow_f + mean_coeff$btempgrowpptgrowsex_f
M_params$flow_pptgrow2<-mean_coeff$bpptgrow2_f + mean_coeff$bpptgrow2sex_f
M_params$flow_pptdorm2<-mean_coeff$bpptdorm2_f + mean_coeff$bpptdorm2sex_f
M_params$flow_tempgrow2<-mean_coeff$btempgrow2_f + mean_coeff$btempgrow2sex_f
M_params$flow_tempdorm2<-mean_coeff$btempdorm2_f + mean_coeff$btempdorm2sex_f
## panicles
F_params$panic_mu <- mean_coefp$b0_p
F_params$panic_size <- mean_coefp$bsize_p
F_params$panic_pptgrow <- mean_coefp$bpptgrow_p 
F_params$panic_pptdorm <- mean_coefp$bpptdorm_p 
F_params$panic_tempgrow <- mean_coefp$btempgrow_p 
F_params$panic_tempdorm <- mean_coefp$btempdorm_p  
F_params$panic_tempdorm_pptdorm<-mean_coefp$btempdormpptdorm_p
F_params$panic_tempgrow_pptgrow<-mean_coefp$btempgrowpptgrow_p
F_params$panic_pptgrow2<-mean_coefp$bpptgrow2_p
F_params$panic_pptdorm2<-mean_coefp$bpptdorm2_p
F_params$panic_tempgrow2<-mean_coefp$btempgrow2_p
F_params$panic_tempdorm2<-mean_coefp$btempdorm2_p
M_params$panic_mu <- mean_coefp$b0_p + mean_coefp$bsex_p  
M_params$panic_size <- mean_coefp$bsize_p + mean_coefp$bsizesex_p
M_params$panic_pptgrow <- mean_coefp$bpptgrow_p  + mean_coefp$bpptgrowsex_p 
M_params$panic_pptdorm <- mean_coefp$bpptdorm_p  + mean_coefp$bpptdormsex_p 
M_params$panic_tempgrow <- mean_coefp$btempgrow_p + mean_coefp$btempgrowsex_p #
M_params$panic_tempdorm <-  mean_coefp$btempdorm_p  + mean_coefp$btempdormsex_p  
M_params$panic_tempdorm_pptdorm<-mean_coefp$btempdormpptdorm_p + mean_coefp$btempdormpptdormsex_p
M_params$panic_tempgrow_pptgrow<-mean_coefp$btempgrowpptgrow_p + mean_coefp$btempgrowpptgrowsex_p
M_params$panic_pptgrow2<-mean_coefp$bpptgrow2_p + mean_coefp$bpptgrow2sex_p
M_params$panic_pptdorm2<-mean_coefp$bpptdorm2_p + mean_coefp$bpptdorm2sex_p
M_params$panic_tempgrow2<-mean_coefp$btempgrow2_p + mean_coefp$btempgrow2sex_p
M_params$panic_tempdorm2<-mean_coefp$btempdorm2_p + mean_coefp$btempdorm2sex_p
## seed viability and misc fertility params
F_params$v0 <- mean_coefv$v0 
F_params$a_v <- mean_coefv$a_v 
F_params$ov_per_inf <- mean_coefv$lambda_d 
F_params$germ <- mean_coefv$m 
F_params$PSR <- 0.5
## use POAU seedling survival for females and males
F_params$sdlg_surv <- M_params$sdlg_surv <- sdlg_surv$sdlg_surv
## set max size equal between the sexes
F_params$max_size <- M_params$max_size <- round(quantile(na.omit((poar.clim_seasonal$tillerN_t1)),probs=0.95))

# Seed viability figure ---- 
viab_dat   <- viabVr %>% 
  select( plot, totS, yesMaybe, sr_f ) %>% 
  rename( SR        = sr_f,
          y_viab = yesMaybe,
          tot_seeds_viab = totS) %>% 
  select(y_viab, tot_seeds_viab, SR ) %>% 
  na.omit

viab_pars <- rstan::extract(fit_full, pars = quote_bare(v0,a_v,m,lambda_d))
post_draws <- sample.int(length(surv_coef$b0_s), n_post_draws)
#print figure
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/seed_viabitility.pdf",useDingbats = F)
plot(jitter(viab_dat$SR,75),jitter((viab_dat$y_viab / viab_dat$tot_seeds_viab),75),
     type="n",xlab="Operational sex ratio",
     ylab="Seed viability",cex.lab=1.2)
for(p in 1:n_post_draws){
  lines(seq(0,1,0.01),
        viab_pars$v0[post_draws[p]] * (1 - seq(0,1,0.01) ^ viab_pars$a_v[post_draws[p]]),
        col=alpha("#E6AB02",0.1))
}
points(jitter(viab_dat$SR,75),jitter((viab_dat$y_viab / viab_dat$tot_seeds_viab),75),
       cex = 5 * (viab_dat$tot_seeds_viab / max(viab_dat$tot_seeds_viab)),lwd=2,col=alpha("#E6AB02",1),bg="#E6AB02",pch=21)
dev.off()



# Sex ratio observed data----
### Climate values across species range
ppt<-terra::rast(clim_1990_2019[[3]])
study_area<-terra::vect("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/USA_vector_polygon/States_shapefile.shp")
study_area <- study_area[(study_area$State_Name %in% c("TEXAS","OKLAHOMA","KANSAS")), ]
plot(study_area)
crop_study_area <- terra::crop(ppt, study_area,mask=TRUE)
crop_study_area_raster<-raster::raster(crop_study_area)
coordinates_study_area<-read.csv("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/coordinate10_study_area.csv")
coordinates(coordinates_study_area) <- ~ x + y
CRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(coordinates_study_area) <- CRS1
plot(crop_study_area_raster)
plot(study_area,add=T)
plot(coordinates_study_area,add=T)

# Past
clim_past_values<-terra::extract(clim_1901_1930,coordinates_study_area)
clim_past_values<-as.data.frame(clim_past_values)
clim_past_values$zpptgrow<-(clim_past_values$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
clim_past_values$ztempgrow<-(clim_past_values$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
clim_past_values$zpptdorm<-(clim_past_values$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
clim_past_values$ztempdorm<-(clim_past_values$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)
clim_past_values<-data.frame(clim_past_values,coordinates_study_area)
dim(clim_past_values)

# Present
clim_current_values<-terra::extract(clim_1990_2019,coordinates_study_area)
clim_current_values<-as.data.frame(clim_current_values)
clim_current_values$zpptgrow<-(clim_current_values$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
clim_current_values$ztempgrow<-(clim_current_values$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
clim_current_values$zpptdorm<-(clim_current_values$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
clim_current_values$ztempdorm<-(clim_current_values$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)
clim_current_values<-data.frame(clim_current_values,coordinates_study_area)
#dim(clim_current_values)

# Future
clim_miroc_85_values<-terra::extract(clim_miroc_85,coordinates_study_area)
climate_miroc85_values<-data.frame(clim_miroc_85_values)
climate_miroc85_values$zpptgrow<-(climate_miroc85_values$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
climate_miroc85_values$ztempgrow<-(climate_miroc85_values$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
climate_miroc85_values$zpptdorm<-(climate_miroc85_values$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
climate_miroc85_values$ztempdorm<-(climate_miroc85_values$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)
climate_miroc85_values<-data.frame(climate_miroc85_values,coordinates_study_area)
# dim(climate_miroc85_values)

clim_miroc_45_values<-terra::extract(clim_miroc_45,coordinates_study_area)
climate_miroc45_values<-data.frame(clim_miroc_45_values)
climate_miroc45_values$zpptgrow<-(climate_miroc45_values$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
climate_miroc45_values$ztempgrow<-(climate_miroc45_values$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
climate_miroc45_values$zpptdorm<-(climate_miroc45_values$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
climate_miroc45_values$ztempdorm<-(climate_miroc45_values$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)
climate_miroc45_values<-data.frame(climate_miroc45_values,coordinates_study_area)
# dim(clim_miroc_45_values)

clim_acc_85_values<-terra::extract(clim_access_85,coordinates_study_area)
climate_acc_85_values<-data.frame(clim_acc_85_values)
climate_acc_85_values$zpptgrow<-(climate_acc_85_values$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
climate_acc_85_values$ztempgrow<-(climate_acc_85_values$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
climate_acc_85_values$zpptdorm<-(climate_acc_85_values$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
climate_acc_85_values$ztempdorm<-(climate_acc_85_values$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)
climate_acc_85_values<-data.frame(climate_acc_85_values,coordinates_study_area)

clim_acc_45_values<-terra::extract(clim_access_45,coordinates_study_area)
climate_acc_45_values<-data.frame(clim_acc_45_values)
climate_acc_45_values$zpptgrow<-(climate_acc_45_values$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
climate_acc_45_values$ztempgrow<-(climate_acc_45_values$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
climate_acc_45_values$zpptdorm<-(climate_acc_45_values$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
climate_acc_45_values$ztempdorm<-(climate_acc_45_values$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)
climate_acc_45_values<-data.frame(climate_acc_45_values,coordinates_study_area)

clim_cmc_85_values<-terra::extract(clim_cmcc_85,coordinates_study_area)
climate_cmc_85_values<-data.frame(clim_cmc_85_values)
climate_cmc_85_values$zpptgrow<-(climate_cmc_85_values$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
climate_cmc_85_values$ztempgrow<-(climate_cmc_85_values$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
climate_cmc_85_values$zpptdorm<-(climate_cmc_85_values$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
climate_cmc_85_values$ztempdorm<-(climate_cmc_85_values$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)
climate_cmc_85_values<-data.frame(climate_cmc_85_values,coordinates_study_area)

clim_cmc_45_values<-terra::extract(clim_cmcc_45,coordinates_study_area)
climate_cmc_45_values<-data.frame(clim_cmc_45_values)
climate_cmc_45_values$zpptgrow<-(climate_cmc_45_values$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
climate_cmc_45_values$ztempgrow<-(climate_cmc_45_values$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
climate_cmc_45_values$zpptdorm<-(climate_cmc_45_values$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
climate_cmc_45_values$ztempdorm<-(climate_cmc_45_values$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)
climate_cmc_45_values<-data.frame(climate_cmc_45_values,coordinates_study_area)

clim_ces_85_values<-terra::extract(clim_cesm_85,coordinates_study_area)
climate_ces_85_values<-data.frame(clim_ces_85_values)
climate_ces_85_values$zpptgrow<-(climate_ces_85_values$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
climate_ces_85_values$ztempgrow<-(climate_ces_85_values$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
climate_ces_85_values$zpptdorm<-(climate_ces_85_values$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
climate_ces_85_values$ztempdorm<-(climate_ces_85_values$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)
climate_ces_85_values<-data.frame(climate_ces_85_values,coordinates_study_area)

clim_ces_45_values<-terra::extract(clim_cesm_45,coordinates_study_area)
climate_ces_45_values<-data.frame(clim_ces_45_values)
climate_ces_45_values$zpptgrow<-(climate_ces_45_values$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
climate_ces_45_values$ztempgrow<-(climate_ces_45_values$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
climate_ces_45_values$zpptdorm<-(climate_ces_45_values$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
climate_ces_45_values$ztempdorm<-(climate_ces_45_values$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)
climate_ces_45_values<-data.frame(climate_ces_45_values,coordinates_study_area)

### Population viability values 
lambbacurrent <- read.csv("https://www.dropbox.com/scl/fi/5xl6vxw456m942qbtynqi/poar_current_1-400K.csv?rlkey=y76suoez6xbk21xkh52iaym9k&dl=1", stringsAsFactors = F) 
lambbapast <- read.csv("https://www.dropbox.com/scl/fi/l1509slfkm5lpecbpscfp/poar_past_1-400K.csv?rlkey=83pxzva68d98qk090kk1kmfkt&dl=1", stringsAsFactors = F)
lambbamiroc45 <- read.csv("https://www.dropbox.com/scl/fi/kq7p7zpnr4n910ed9lf1n/poar_miroc45_1-400K.csv?rlkey=0381obcnwjhzjkcx5kj8krmua&dl=1", stringsAsFactors = F)
lambbamiroc85 <- read.csv("https://www.dropbox.com/scl/fi/tzat8n3x9dgj4mr2kj2a0/poar_miroc85_1-400K.csv?rlkey=yv10l7v2ii4ljossoviebyp59&dl=1", stringsAsFactors = F)
lambbacmc45 <- read.csv("https://www.dropbox.com/scl/fi/ekjo0ezv0gn2pfoj69blt/poar_cmc45_1-400K.csv?rlkey=25x99ojbmc021qswxl2ehzq3t&dl=1", stringsAsFactors = F)
lambbacmc85 <- read.csv("https://www.dropbox.com/scl/fi/299bagx6imcgm534bcmhs/poar_cmc85_1-400K.csv?rlkey=fd4qs8grnkyua7o99w0oz9pyt&dl=1", stringsAsFactors = F)
lambbaacc45 <- read.csv("https://www.dropbox.com/scl/fi/973dwzsk6y7hzpoov5ljo/poar_acc45_1-400K.csv?rlkey=qsib01gvqjjgmye7pofnmvm7n&dl=1", stringsAsFactors = F)
lambbaacc85 <- read.csv("https://www.dropbox.com/scl/fi/ybt7vyi1693nselqqkxhc/poar_acc85_1-400K.csv?rlkey=pfzooz6p1p0bw2dnvydm2gdzf&dl=1", stringsAsFactors = F)
lambbaces45 <- read.csv("https://www.dropbox.com/scl/fi/8ml4oep2lssgg4xd3eon6/poar_ces45_1-400K.csv?rlkey=71ow2cfxcs6uk7iyqfcdijn1y&dl=1", stringsAsFactors = F)
lambbaces85 <- read.csv("https://www.dropbox.com/scl/fi/84d3crazogxfj5p2id4to/poar_ces85_1-400K.csv?rlkey=raur7uo6w3flwippbrom69uqb&dl=1", stringsAsFactors = F)

### Female dominant model
geo_lambbacurrent_fd <- read.csv("https://www.dropbox.com/scl/fi/7zmsid3dxszv3l3tyvoox/lambda_post_current.csv?rlkey=hwy253ls9imsupkzqrduwp7r8&dl=1", stringsAsFactors = F) 
geo_lambbapast_fd <- read.csv("https://www.dropbox.com/scl/fi/9x9gbdv1vefdyd292726c/data_past_fd.csv?rlkey=besdw6mj9ab1nubx18if54w9z&dl=1", stringsAsFactors = F)
geo_lambda_miroc45_fd<-read.csv("https://www.dropbox.com/scl/fi/dmfb2mdhrv5k8otcrc33m/lambda_miroc45.csv?rlkey=1d4j1blbmlx9i3rp61c35g7f8&dl=1", stringsAsFactors = F)
geo_lambba_miroc85_fd <- read.csv("https://www.dropbox.com/scl/fi/a1rmqszey6y9s88887fqq/datamiroc85fd.csv?rlkey=1eckoajb6s4ee6dq61f2o5ymj&dl=1", stringsAsFactors = F)
geo_lambba_cmc45_fd <- read.csv("https://www.dropbox.com/scl/fi/oarso87jb6q5c30gsqm3k/lambda_cmcc45_fd.csv?rlkey=w7xlibvxnc2etc3iliog960dq&dl=1", stringsAsFactors = F)
geolambba_cmc85_fd <- read.csv("https://www.dropbox.com/scl/fi/d8uv2mbb2g47ywly0cz63/data_cmc85_fd.csv?rlkey=o4wcrp7hkqmdyow19eue8e1dh&dl=1", stringsAsFactors = F)
geo_lambba_acc45_fd <- read.csv("https://www.dropbox.com/scl/fi/p7xkr30t6zf24gtn4adbl/data_acc45_fd.csv?rlkey=2qv22qe7p3bahmn24mmtff1rm&dl=1", stringsAsFactors = F)
geolambba_acc85_fd <- read.csv("https://www.dropbox.com/scl/fi/5ged79dtckp3z2k570sm8/data_acc85_fd.csv?rlkey=l97zn5wfci6icbbvv8gcj6z5c&dl=1", stringsAsFactors = F)
geo_lambba_ces45_fd <- read.csv("https://www.dropbox.com/scl/fi/zdj4ra2qaghb4y16uiuow/data_ces45_fd.csv?rlkey=wl8t4ih9s48dr674xti0wwwww&dl=1", stringsAsFactors = F)
geolambba_ces85_fd <- read.csv("https://www.dropbox.com/scl/fi/e4mksgd5hpayfwd3oy0sz/data_ces85_fd.csv?rlkey=miacot7xh5yan0vk6nbwk9q1o&dl=1", stringsAsFactors = F)

## Estimating Probability of lambda being higher than 1 across the landscape---- 
### Past conditions-----
lambbapast %>% 
  unique %>% 
  pivot_wider(
    names_from = posterior_id, 
    values_from = lambda_diff)->pivot_lambda_past
#dim(pivot_lambda_past) 

pivot_lambda_past_final<-pivot_lambda_past[,-c(1:8)]
prob_lambda_past_post<-c() 
for(l in 1:nrow(pivot_lambda_past_final)){
  prob_lambda_past_post[l]<-mean(pivot_lambda_past_final[l,]>1,na.rm=T)
}
lam_prob_past<-data.frame(pivot_lambda_past[,1:8],Prlambda=prob_lambda_past_post)
geo_prlambda_past <- merge(x = lam_prob_past,y =clim_past_values,by=c("pptdorm","tempdorm"), all.x = TRUE) 

prob_lambda_past_fd<-c()
for(l in 1:ncol(geo_lambbapast_fd)){
  prob_lambda_past_fd[l]<-mean(geo_lambbapast_fd[,l]>1,na.rm=T)
}

clim_past_values_clean<-na.omit(clim_past_values)
lam_prob_past_fd<-data.frame(clim_past_values_clean[,9:10],Prlambda=prob_lambda_past_fd)

### Present conditions
lambbacurrent %>% 
  unique %>% 
  pivot_wider( values_from = 'lambda_diff',
               names_from  = 'posterior_id' )->pivot_lambda_current
dim(pivot_lambda_current)
pivot_lambda_current_final<-pivot_lambda_current[,-c(1:8)]
# dim(pivot_lambda_current_final)
# Calculate probability of lambda being higher 1
prob_lambda_current_post<-c()
for(l in 1:nrow(pivot_lambda_current_final)){
  prob_lambda_current_post[l]<-mean(pivot_lambda_current_final[l,]>1,na.rm=T)
}
lam_prob_current<-data.frame(pivot_lambda_current[,1:8],Prlambda=prob_lambda_current_post)
geo_prlambda_current <- merge(x =lam_prob_current ,y = clim_current_values,by=c("pptdorm","tempdorm"),all.x=TRUE) 

prob_lambda_current_fd<-c()
for(l in 1:ncol(geo_lambbacurrent_fd)){
  prob_lambda_current_fd[l]<-mean(geo_lambbacurrent_fd[,l]>1,na.rm=T)
}
lam_prob_current_fd<-data.frame(clim_current_values[,9:10],Prlambda=prob_lambda_current_fd)

### Future conditions----
#### MIROC 5----
lambbamiroc45 %>% 
  unique %>% 
  pivot_wider( values_from = 'lambda_diff',
               names_from  = 'posterior_id' )->pivot_lambda_miroc45
# dim(pivot_lambda_miroc45)
pivot_lambda_miroc45_final<-pivot_lambda_miroc45[,-c(1:8)]
# dim(pivot_lambda_45_final)
prob_lambda_miroc45_post<-c()
for(l in 1:nrow(pivot_lambda_miroc45_final)){
  prob_lambda_miroc45_post[l]<-mean(pivot_lambda_miroc45_final[l,]>1,na.rm=T)
}
lam_prob_miroc45<-data.frame(pivot_lambda_miroc45[,1:8],Prlambda=prob_lambda_miroc45_post)
geo_prlambda_miroc45 <- merge(x = lam_prob_miroc45,y =climate_miroc45_values,by=c("pptdorm","tempdorm"),all.x=TRUE) 

lambbamiroc85 %>% 
  unique %>% 
  pivot_wider( values_from = 'lambda_diff',
               names_from  = 'posterior_id' )->pivot_lambda_miroc85
# dim(pivot_lambda_miroc85)
pivot_lambda_miroc85_final<-pivot_lambda_miroc85[,-c(1:8)]
# dim(pivot_lambda_miroc85_final)
prob_lambda_miroc85_post<-c()
for(l in 1:nrow(pivot_lambda_miroc85_final)){
  prob_lambda_miroc85_post[l]<-mean(pivot_lambda_miroc85_final[l,]>1,na.rm=T)
}
lam_prob_miroc85<-data.frame(pivot_lambda_miroc85[,1:8],Prlambda=prob_lambda_miroc85_post)
geo_prlambda_miroc85 <- merge(x = lam_prob_miroc85,y =climate_miroc85_values,by=c("pptdorm","tempdorm"))

prob_lambda_miroc45_fd<-c()
for(l in 1:ncol(geo_lambda_miroc45_fd)){
  prob_lambda_miroc45_fd[l]<-mean(geo_lambda_miroc45_fd[,l]>1,na.rm=T)
}
lam_prob_miroc45_fd<-data.frame(climate_miroc45_values[,9:10],Prlambda=prob_lambda_miroc45_fd)
prob_lambda_miroc85_fd<-c()
dim(geo_lambba_miroc85_fd)
for(l in 1:ncol(geo_lambba_miroc85_fd)){
  prob_lambda_miroc85_fd[l]<-mean(geo_lambba_miroc85_fd[,l]>1,na.rm=T)
}
lam_prob_miroc85_fd<-data.frame(climate_miroc85_values[,9:10],Prlambda=prob_lambda_miroc85_fd)

#### ACCESS ----
lambbaacc45 %>% 
  unique %>% 
  pivot_wider( values_from = 'lambda_diff',
               names_from  = 'posterior_id' )->pivot_lambdaacc_45
# dim(pivot_lambdaacc_45)
pivot_lambdaacc_45_final<-pivot_lambdaacc_45[,-c(1:8)]
prob_lambdaacc_45_post<-c()
for(l in 1:nrow(pivot_lambdaacc_45_final)){
  prob_lambdaacc_45_post[l]<-mean(pivot_lambdaacc_45_final[l,]>1,na.rm=T)
}
lambaacc_prob_45<-data.frame(pivot_lambdaacc_45[,1:8],Prlambda=prob_lambdaacc_45_post)
geo_prlambdaacc_45 <- merge(x = lambaacc_prob_45,y =climate_acc_45_values,by=c("pptdorm","tempdorm"),all.x=TRUE) 

lambbaacc85 %>% 
  unique %>% 
  pivot_wider( values_from = 'lambda_diff',
               names_from  = 'posterior_id' )->pivot_lambdacc_85
pivot_lambdacc_85_final<-pivot_lambdacc_85[,-c(1:8)]
prob_lambdaacc_85_post<-c()
for(l in 1:nrow(pivot_lambdacc_85_final)){
  prob_lambdaacc_85_post[l]<-mean(pivot_lambdacc_85_final[l,]>1,na.rm=T)
}
lambdaacc_prob_85<-data.frame(pivot_lambdacc_85[,1:8],Prlambda=prob_lambdaacc_85_post)
geo_prlambdaacc_85 <- merge(x = lambdaacc_prob_85,y =climate_acc_85_values,by=c("pptdorm","tempdorm")) 

prob_lambda_acc45_fd<-c()
for(l in 1:ncol(geo_lambba_acc45_fd)){
  prob_lambda_acc45_fd[l]<-mean(geo_lambba_acc45_fd[,l]>1,na.rm=T)
}
lam_prob_acc45_fd<-data.frame(climate_acc_45_values[,9:10],Prlambda=prob_lambda_acc45_fd)
prob_lambda_acc85_fd<-c()
dim(geolambba_acc85_fd)
for(l in 1:ncol(geolambba_acc85_fd)){
  prob_lambda_acc85_fd[l]<-mean(geolambba_acc85_fd[,l]>1,na.rm=T)
}
lam_prob_acc85_fd<-data.frame(climate_acc_85_values[,9:10],Prlambda=prob_lambda_acc85_fd)

#### CMC----
lambbacmc45 %>% 
  unique %>% 
  pivot_wider( values_from = 'lambda_diff',
               names_from  = 'posterior_id' )->pivot_lambdacmc_45
pivot_lambdacmc_45_final<-pivot_lambdacmc_45[,-c(1:8)]
prob_lambdacmc_45_post<-c()
for(l in 1:nrow(pivot_lambdacmc_45_final)){
  prob_lambdacmc_45_post[l]<-mean(pivot_lambdacmc_45_final[l,]>1,na.rm=T)
}
lamcmc_prob_45<-data.frame(pivot_lambdacmc_45[,1:8],Prlambda=prob_lambdacmc_45_post)

geo_prlambdacmc_45 <- merge(x = lamcmc_prob_45,y =climate_cmc_45_values,by=c("pptdorm","tempdorm"),all.x=TRUE) 

lambbacmc85 %>% 
  unique %>% 
  pivot_wider( values_from = 'lambda_diff',
               names_from  = 'posterior_id' )->pivot_lambdacmc_85
# dim(pivot_lambda_current)
pivot_lambdacmc_85_final<-pivot_lambdacmc_85[,-c(1:8)]
dim(pivot_lambdacmc_85_final)
# Calculate probability of lambda above 1
prob_lambdacmc_85_post<-c()
for(l in 1:nrow(pivot_lambdacmc_85_final)){
  prob_lambdacmc_85_post[l]<-mean(pivot_lambdacmc_85_final[l,]>1,na.rm=T)
}
lamcmc_prob_85<-data.frame(pivot_lambdacmc_85[,1:8],Prlambda=prob_lambdacmc_85_post)
geo_prlambdacmc_85 <- merge(x = lamcmc_prob_85,y =climate_cmc_85_values,by=c("pptdorm","tempdorm")) # merge the demographic data with climatic data for each
prob_lambda_cmc45_fd<-c()
dim(geo_lambba_cmc45_fd)
for(l in 1:ncol(geo_lambba_cmc45_fd)){
  prob_lambda_cmc45_fd[l]<-mean(geo_lambba_cmc45_fd[,l]>1,na.rm=T)
}
lam_prob_cmc45_fd<-data.frame(climate_cmc_45_values[,9:10],Prlambda=prob_lambda_cmc45_fd)
prob_lambda_cmc85_fd<-c()
for(l in 1:ncol(geolambba_cmc85_fd)){
  prob_lambda_cmc85_fd[l]<-mean(geolambba_cmc85_fd[,l]>1,na.rm=T)
}
lam_prob_cmc85_fd<-data.frame(climate_cmc_85_values[,9:10],Prlambda=prob_lambda_cmc85_fd)

#### CESM----
lambbaces45 %>% 
  unique %>% 
  pivot_wider( values_from = 'lambda_diff',
               names_from  = 'posterior_id' )->pivot_lambdaces_45
pivot_lambdaces_45_final<-pivot_lambdaces_45[,-c(1:8)]
prob_lambdaces_45_post<-c()
for(l in 1:nrow(pivot_lambdaces_45_final)){
  prob_lambdaces_45_post[l]<-mean(pivot_lambdaces_45_final[l,]>1,na.rm=T)
}
lambaces_prob_45<-data.frame(pivot_lambdaces_45[,1:8],Prlambda=prob_lambdaces_45_post)
geo_prlambdaces_45 <- merge(x = lambaces_prob_45,y =climate_ces_45_values,by=c("pptdorm","tempdorm"),all.x=TRUE)
lambbaces85 %>% 
  unique %>% 
  pivot_wider( values_from = 'lambda_diff',
               names_from  = 'posterior_id' )->pivot_lambdaces_85
pivot_lambdaces_85_final<-pivot_lambdaces_85[,-c(1:8)]
prob_lambdaces_85_post<-c()
for(l in 1:nrow(pivot_lambdaces_85_final)){
  prob_lambdaces_85_post[l]<-mean(pivot_lambdaces_85_final[l,]>1,na.rm=T)
}
lambdaces_prob_85<-data.frame(pivot_lambdaces_85[,1:8],Prlambda=prob_lambdaces_85_post)
geo_prlambdaces_85 <- merge(x = lambdaces_prob_85,y =climate_ces_85_values,by=c("pptdorm","tempdorm")) # merge the demographic data with climatic data for each
prob_lambda_ces45_fd<-c()
dim(geo_lambba_ces45_fd)
for(l in 1:ncol(geo_lambba_ces45_fd)){
  prob_lambda_ces45_fd[l]<-mean(geo_lambba_ces45_fd[,l]>1,na.rm=T)
}
lam_prob_ces45_fd<-data.frame(climate_ces_45_values[,9:10],Prlambda=prob_lambda_ces45_fd)
prob_lambda_ces85_fd<-c()
dim(geolambba_ces85_fd)
for(l in 1:ncol(geolambba_ces85_fd)){
  prob_lambda_ces85_fd[l]<-mean(geolambba_ces85_fd[,l]>1,na.rm=T)
}
lam_prob_ces85_fd<-data.frame(climate_ces_85_values[,9:10],Prlambda=prob_lambda_ces85_fd)

# Sex ratio derived from models----
# Importing sex ratio
# Current
zclim_current = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_current)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_current$zpptgrow<-(clim_current$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_current$ztempgrow<-(clim_current$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_current$zpptdorm<-(clim_current$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_current$ztempdorm<-(clim_current$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

## Past
zclim_past = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_past)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_past$zpptgrow<-(clim_past$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_past$ztempgrow<-(clim_past$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_past$zpptdorm<-(clim_past$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_past$ztempdorm<-(clim_past$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

## MIROC45
zclim_miroc45 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_miroc45)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_miroc45$zpptgrow<-(clim_miroc45$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_miroc45$ztempgrow<-(clim_miroc45$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_miroc45$zpptdorm<-(clim_miroc45$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_miroc45$ztempdorm<-(clim_miroc45$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

## ACCESS45
zclim_acc45 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_acc45)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_acc45$zpptgrow<-(clim_access45$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_acc45$ztempgrow<-(clim_access45$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_acc45$zpptdorm<-(clim_access45$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_acc45$ztempdorm<-(clim_access45$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

## CMC45
zclim_cmc45 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_cmc45)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_cmc45$zpptgrow<-(clim_cmcc45$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_cmc45$ztempgrow<-(clim_cmcc45$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_cmc45$zpptdorm<-(clim_cmcc45$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_cmc45$ztempdorm<-(clim_cmcc45$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

## CES45
zclim_ces45 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_ces45)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_ces45$zpptgrow<-(clim_cesm45$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_ces45$ztempgrow<-(clim_cesm45$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_ces45$zpptdorm<-(clim_cesm45$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_ces45$ztempdorm<-(clim_cesm45$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

## Gather the data in to a single dataset
data_site_pptgrow45<-data.frame(zclim_miroc45$zpptgrow,zclim_acc45$zpptgrow,zclim_cmc45$zpptgrow,zclim_ces45$zpptgrow)
data_site_tempgrow45<-data.frame(zclim_miroc45$ztempgrow,zclim_acc45$ztempgrow,zclim_cmc45$ztempgrow,zclim_ces45$ztempgrow)
data_site_pptdorm45<-data.frame(zclim_miroc45$zpptdorm,zclim_acc45$zpptdorm,zclim_cmc45$zpptdorm,zclim_ces45$zpptdorm)
data_site_tempdorm45<-data.frame(zclim_miroc45$ztempdorm,zclim_acc45$ztempdorm,zclim_cmc45$ztempdorm,zclim_ces45$ztempdorm)
data_site_45<-data.frame(zpptgrow=rowMeans(data_site_pptgrow45),ztempgrow=rowMeans(data_site_tempgrow45),zpptdorm=rowMeans(data_site_pptdorm45),ztempdorm=rowMeans(data_site_tempdorm45))

## MIROC85

zclim_miroc85 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_miroc85)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_miroc85$zpptgrow<-(clim_miroc85$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_miroc85$ztempgrow<-(clim_miroc85$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_miroc85$zpptdorm<-(clim_miroc85$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_miroc85$ztempdorm<-(clim_miroc85$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

## ACCESS85
zclim_acc85 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_acc85)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_acc85$zpptgrow<-(clim_access85$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_acc85$ztempgrow<-(clim_access85$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_acc85$zpptdorm<-(clim_access85$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_acc85$ztempdorm<-(clim_access85$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

## CMC85
zclim_cmc85 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_cmc85)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_cmc85$zpptgrow<-(clim_cmcc85$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_cmc85$ztempgrow<-(clim_cmcc85$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_cmc85$zpptdorm<-(clim_cmcc85$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_cmc85$ztempdorm<-(clim_cmcc85$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

## CES85
zclim_ces85 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_ces85)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_ces85$zpptgrow<-(clim_cesm85$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_ces85$ztempgrow<-(clim_cesm85$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_ces85$zpptdorm<-(clim_cesm85$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_ces85$ztempdorm<-(clim_cesm85$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

data_site_pptgrow85<-data.frame(zclim_miroc85$zpptgrow,zclim_acc85$zpptgrow,zclim_cmc85$zpptgrow,zclim_ces85$zpptgrow)
data_site_tempgrow85<-data.frame(zclim_miroc85$ztempgrow,zclim_acc85$ztempgrow,zclim_cmc85$ztempgrow,zclim_ces85$ztempgrow)
data_site_pptdorm85<-data.frame(zclim_miroc85$zpptdorm,zclim_acc85$zpptdorm,zclim_cmc85$zpptdorm,zclim_ces85$zpptdorm)
data_site_tempdorm85<-data.frame(zclim_miroc85$ztempdorm,zclim_acc85$ztempdorm,zclim_cmc85$ztempdorm,zclim_ces85$ztempdorm)
data_site_85<-data.frame(zpptgrow=rowMeans(data_site_pptgrow85),ztempgrow=rowMeans(data_site_tempgrow85),zpptdorm=rowMeans(data_site_pptdorm85),ztempdorm=rowMeans(data_site_tempdorm85))


data_past_present_45_85<-rbind(zclim_past,zclim_current,data_site_45,data_site_85)


sexratio_dormant<-read_csv(url("https://www.dropbox.com/scl/fi/i4e5amud3v0td6jz10j86/SR_OSR_dorm.csv?rlkey=bk1gyfcjn68umysxazydgqrmo&dl=1"))
all_pptdorm_seq<-seq(min(data_past_present_45_85$zpptdorm*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm)),max(data_past_present_45_85$zpptdorm*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm)),length.out=30)
all_tempdorm_seq<-seq(min(data_past_present_45_85$ztempdorm*sd(poar_2015_2016$tempdorm)+mean(poar_2015_2016$tempdorm)),max(data_past_present_45_85$ztempdorm*sd(poar_2015_2016$tempdorm)+mean(poar_2015_2016$tempdorm)),length.out=30)

OSR_mtrx_dorm <- matrix(sexratio_dormant$OSR, nrow = 30, dimnames = list(all_pptdorm_seq,all_tempdorm_seq))

sexratio_growing <- read_csv(url("https://www.dropbox.com/scl/fi/vr2pkz8ot2v6azqf18y0t/OSR_grow.csv?rlkey=l1h2re7kf7q155fuwg8m5m8ot&dl=1"))
all_pptgrow_seq<-seq(min(data_past_present_45_85$zpptgrow*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow)),max(data_past_present_45_85$zpptgrow*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow)),length.out=30)
all_tempgrow_seq<-seq(min(data_past_present_45_85$ztempgrow*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow)),max(data_past_present_45_85$ztempgrow*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow)),length.out=30)
OSR_mtrx_grow <- matrix(sexratio_growing$OSR, nrow = 30, dimnames = list(all_pptgrow_seq,all_tempgrow_seq))

scal_breaks_dorm <- c(seq(0, 1, length.out = 101))
scal_breaks_grow <- c(seq(0.5, 1, length.out = 101))
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/OSR.pdf",width=5,height=7,useDingbats = F)
par(mar=c(5,5,1,1),mfrow=c(2,1))
fields::image.plot(all_pptdorm_seq,all_tempdorm_seq,na.omit(OSR_mtrx_dorm),col=topo.colors(100),xlab="Precipitation dormant",ylab="Temperature dormant",main="",breaks = scal_breaks) 
contour(all_pptdorm_seq,all_tempdorm_seq,OSR_mtrx_dorm,add=T,labcex=0.75,col="black") 
mtext("A",side = 3, adj = 0,cex=1.25)
text(x=x_current_d,y=y_current_d,"+",col="black",cex=1.5)
text(x=x_past_d,y=y_past_d,"o",col="black")
text(x=x_miroc45_d,y=y_miroc45_d,"*",col="black",cex=2)
text(x=x_miroc85_d,y=y_miroc85_d,".",col="black",cex=4.5)
fields::image.plot(all_pptgrow_seq,all_tempgrow_seq,OSR_mtrx_grow,col=topo.colors(100),xlab="Precipitation growing",ylab="Temperature growing",main="")
contour(all_pptgrow_seq,all_tempgrow_seq,OSR_mtrx_grow,add=T,labcex=0.75,col="black")
mtext("B",side = 3, adj = 0,cex=1.25)
# title(main="B" ,adj=0,cex.main=1.4,font=1)
text(x=x_current_g,y=y_current_g,"+",col="black",cex=1.5)
text(x=x_past_g,y=y_past_g,"o",col="black")
text(x=x_miroc45_g,y=y_miroc45_g,"*",col="black",cex=2)
text(x=x_miroc85_g,y=y_miroc85_g,".",col="black",cex=4.5)
text(x=x_miroc85_g,y=y_miroc85_g,".",col="black",cex=4.5)
dev.off()

x_current_d<-clim_current$pptdorm
y_current_d<-clim_current$tempdorm
x_past_d<-clim_past$pptdorm
y_past_d<-clim_past$tempdorm
x_miroc45_d<-clim_miroc45$pptdorm
y_miroc45_d<-clim_miroc45$tempdorm
x_miroc85_d<-clim_miroc85$pptdorm
y_miroc85_d<-clim_miroc85$tempdorm

x_current_g<-clim_current$pptgrow
y_current_g<-clim_current$tempgrow
x_past_g<-clim_past$pptgrow
y_past_g<-clim_past$tempgrow
x_miroc45_g<-clim_miroc45$pptgrow
y_miroc45_g<-clim_miroc45$tempgrow
x_miroc85_g<-clim_miroc85$pptgrow
y_miroc85_g<-clim_miroc85$tempgrow



pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/OSR.pdf",width=5,height=7,useDingbats = F)
par(mar=c(5,5,1,1),mfrow=c(2,1))
fields::image.plot(all_pptdorm_seq,all_tempdorm_seq,na.omit(OSR_mtrx_dorm),col=topo.colors(100),xlab="Dormant season precip",ylab="Dormant season temp",main="",breaks = scal_breaks_dorm,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,legend.mar = 4,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text="OSR", side=3, font=3, line=0.3, cex=0.6)) 
# contour(all_pptdorm_seq,all_tempdorm_seq,OSR_mtrx_dorm,add=T,labcex=0.75,col="black") 
text(x=x_current_d,y=y_current_d,label = paste0("+"),col="black",cex=1)
text(x=x_past_d,y=y_past_d,label = paste0("o"),col="black",cex=1)
text(x=x_miroc45_d,y=y_miroc45_d,label = paste0("*"),col="black",cex=1)
text(x=x_miroc85_d,y=y_miroc85_d,label = paste0("-"),col="black",cex=1)
# mtext("Survival",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "A",side = 3, adj = 0,cex=1.2)
fields::image.plot(all_pptgrow_seq,all_tempgrow_seq,OSR_mtrx_grow,col=topo.colors(100),xlab="Growing season precip",ylab="Growing season temp",main="",breaks = scal_breaks_grow,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,legend.mar = 4,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text="OSR", side=3, font=3, line=0.3, cex=0.6)) 
# contour(all_pptgrow_seq,all_tempgrow_seq,OSR_mtrx_grow,add=T,labcex=0.75,col="black") 
text(x=x_current_g,y=y_current_g,label = paste0("+"),col="black",cex=1)
text(x=x_past_g,y=y_past_g,label = paste0("o"),col="black",cex=1)
text(x=x_miroc45_g,y=y_miroc45_g,label = paste0("*"),col="black",cex=1)
text(x=x_miroc85_g,y=y_miroc85_g,label = paste0("-"),col="black",cex=1)
# mtext("Survival",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "B",side = 3, adj = 0,cex=1.2)
dev.off()


# Importing the data on Proportion of female (SR) and proportion of female panicles (OSR)
sexratio_past<-read_csv(url("https://www.dropbox.com/scl/fi/yfv9b4w857wg8hcnasg6b/future_spatial_past.csv?rlkey=1q5xovdf33g2woql3p56iq9dr&dl=1"))
sexratio_current<-read_csv(url("https://www.dropbox.com/scl/fi/92rsr3ek3pzmvg4ej3pyp/future_spatial_current.csv?rlkey=xqbw7cei5u6wv2te41bzdohs0&dl=1"))
sexratio_cmc45<-read_csv(url("https://www.dropbox.com/scl/fi/cxalajldmu9dngs8767xv/future_spatial_cmc45.csv?rlkey=zgnl8kpu78bml8pw9fl0x81a7&dl=1"))
sexratio_cmc85<-read_csv(url("https://www.dropbox.com/scl/fi/jhof247ui1eor19iukswi/future_spatial_cmc85.csv?rlkey=w9b7fx3y5jxmvqyg5cq792pim&dl=1"))
sexratio_acc45<-read_csv(url("https://www.dropbox.com/scl/fi/5hhvjv3ikv95vjuaci72m/future_spatial_acc45.csv?rlkey=7tq5wwdy23uwre8n2ied823mj&dl=1"))
sexratio_acc85<-read_csv(url("https://www.dropbox.com/scl/fi/uf9d2c1eoj7a5sd3w2jsa/future_spatial_acc85.csv?rlkey=ntnhaa7dfzbmizqqsgfuslhl0&dl=1"))
sexratio_ces45<-read_csv(url("https://www.dropbox.com/scl/fi/4gqeocp3sgnqup1htzci4/future_spatial_ces45.csv?rlkey=oxpyohr6ghxttc2kxviswwh7x&dl=1"))
sexratio_ces85<-read_csv(url("https://www.dropbox.com/scl/fi/wdqmv4fsg9xdp50quesdz/future_spatial_ces85.csv?rlkey=qda9kjbsjf55r6w6vklf9s099&dl=1"))
sexratio_miroc45<-read_csv(url("https://www.dropbox.com/s/s7ps4029miuefnh/future_spatial_miroc45.csv?dl=1"))
sexratio_miroc85<-read_csv(url("https://www.dropbox.com/scl/fi/y9pf380aw5103h6l81dhi/future_spatial_miroc85.csv?rlkey=a7qk5qj27liyopkn54up4va0b&dl=1"))


## SR and OSR for Pr lambda higher than 1 
sr_past <- merge(x = sexratio_past,y =lam_prob_past,by=c("pptdorm","tempdorm"))
sr_current <- merge(x = sexratio_current,y =lam_prob_current,by=c("pptdorm","tempdorm")) 
sr_cmc_85 <- merge(x = sexratio_cmc85,y =lamcmc_prob_85,by=c("pptdorm","tempdorm"))
sr_cmc_45 <- merge(x = sexratio_cmc45,y =lamcmc_prob_45,by=c("pptdorm","tempdorm")) 
sr_miroc_85 <- merge(x = sexratio_miroc85,y =lam_prob_miroc85,by=c("pptdorm","tempdorm"))
sr_miroc_45 <- merge(x = sexratio_miroc45,y =lam_prob_miroc45,by=c("pptdorm","tempdorm")) 
sr_ces_85 <- merge(x = sexratio_ces85,y =lambdaces_prob_85,by=c("pptdorm","tempdorm"))
sr_ces_45 <- merge(x = sexratio_ces45,y =lambaces_prob_45,by=c("pptdorm","tempdorm")) 
sr_acc_85 <- merge(x =sexratio_acc85 ,y =lambdaacc_prob_85,by=c("pptdorm","tempdorm"))
sr_acc_45 <- merge(x = sexratio_acc45,y = lambaacc_prob_45,by=c("pptdorm","tempdorm")) 

# Define sex ratio as a function of Pr(lambda>1)
sr_past$OSR_final <- ifelse(sr_past$Prlambda >0.5, sr_past$OSR, NA)
sr_past$SR_final <- ifelse(sr_past$Prlambda >0.5, sr_past$SR, NA)

sr_current$OSR_final <- ifelse(sr_current$Prlambda >0.5, sr_current$OSR, NA)
sr_current$SR_final <- ifelse(sr_current$Prlambda >0.5, sr_current$SR, NA)

sr_cmc_85$OSR_final <- ifelse(sr_cmc_85$Prlambda >0.5, sr_cmc_85$OSR, NA)
sr_cmc_85$SR_final <- ifelse(sr_cmc_85$Prlambda >0.5, sr_cmc_85$SR, NA)

sr_cmc_45$OSR_final <- ifelse(sr_cmc_45$Prlambda >0.5, sr_cmc_45$OSR, NA)
sr_cmc_45$SR_final <- ifelse(sr_cmc_45$Prlambda >0.5, sr_cmc_45$SR, NA)

sr_miroc_85$OSR_final <- ifelse(sr_miroc_85$Prlambda >0.5, sr_miroc_85$OSR, NA)
sr_miroc_85$SR_final <- ifelse(sr_miroc_85$Prlambda >0.5, sr_miroc_85$SR, NA)

sr_miroc_45$OSR_final <- ifelse(sr_miroc_45$Prlambda >0.5, sr_miroc_45$OSR, NA)
sr_miroc_45$SR_final <- ifelse(sr_miroc_45$Prlambda >0.5, sr_miroc_45$SR, NA)

sr_ces_85$OSR_final <- ifelse(sr_ces_85$Prlambda >0.5, sr_ces_85$OSR, NA)
sr_ces_85$SR_final <- ifelse(sr_ces_85$Prlambda >0.5, sr_ces_85$SR, NA)

sr_ces_45$OSR_final <- ifelse(sr_ces_45$Prlambda >0.5, sr_ces_45$OSR, NA)
sr_ces_45$SR_final <- ifelse(sr_ces_45$Prlambda >0.5, sr_ces_45$SR, NA)

sr_acc_85$OSR_final <- ifelse(sr_acc_85$Prlambda >0.5, sr_acc_85$OSR, NA)
sr_acc_85$SR_final <- ifelse(sr_acc_85$Prlambda >0.5, sr_acc_85$SR, NA)

sr_acc_45$OSR_final <- ifelse(sr_acc_45$Prlambda >0.5, sr_acc_45$OSR, NA)
sr_acc_45$SR_final <- ifelse(sr_acc_45$Prlambda >0.5, sr_acc_45$SR, NA)

## Density plots----
d_past_OSR <- density(unique_no.NA(sr_past$OSR_final))
d_current_OSR <- density(unique_no.NA(sr_current$OSR_final))
d_cmc45_OSR <- density(unique_no.NA(sr_cmc_45$OSR_final))
d_cmc85_OSR <- density(unique_no.NA(sr_cmc_85$OSR_final))

d_ces45_OSR <- density(unique_no.NA(sr_ces_45$OSR_final))
d_ces85_OSR <- density(unique_no.NA(sr_ces_85$OSR_final))


d_acc45_OSR <- density(unique_no.NA(sr_acc_45$OSR_final))
d_acc85_OSR <- density(unique_no.NA(sr_acc_85$OSR_final))

d_miroc45_OSR <- density(unique_no.NA(sr_miroc_45$OSR_final))
d_miroc85_OSR <- density(unique_no.NA(sr_miroc_85$OSR_final))

d_past_SR <- density(unique_no.NA(sr_past$SR_final))
d_current_SR <- density(unique_no.NA(sr_current$SR_final))
d_cmc45_SR <- density(unique_no.NA(sr_cmc_45$SR_final))
d_cmc85_SR <- density(unique_no.NA(sr_cmc_85$SR_final))

d_ces45_SR <- density(unique_no.NA(sr_ces_45$SR_final))
d_ces85_SR <- density(unique_no.NA(sr_ces_85$SR_final))


d_acc45_SR <- density(unique_no.NA(sr_acc_45$SR_final))
d_acc85_SR <- density(unique_no.NA(sr_acc_85$SR_final))

d_miroc45_SR <- density(unique_no.NA(sr_miroc_45$SR_final))
d_miroc85_SR <- density(unique_no.NA(sr_miroc_85$SR_final))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/POAR_OSR.pdf",width=5,height=5,useDingbats = F)
par(mar=c(5,5,1.5,0.5))
plot(d_past_OSR, lwd = 2, main = "", xlab = "Operational Sex Ratio",xlim=c(0,1),ylim=c(0,6.5),
     col = "#7570B3")
# mtext( "A",side = 3, adj = 0,cex=1.25)
abline(v=mean(sr_past$OSR_final,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_OSR, lwd = 2,col = "#0072B2")
abline(v=mean(sr_current$OSR_final,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_cmc45_OSR, lwd = 2,col = "#E6AB02")
abline(v=mean(sr_cmc_45$OSR_final,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_cmc85_OSR, lwd = 2,col = "#E7298A")
abline(v=mean(sr_cmc_85$OSR_final,na.rm=TRUE),lty=2,col="#E7298A")

legend(0.1,6,bty="n",legend=c("Past","Present","RCP 4.5","RCP 8.5"),lwd=c(2,2,2,2),lty=1,col=c("#7570B3","#0072B2","#E6AB02","#E7298A"),cex=0.8)

dev.off()

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/POAR_OSR_MIROC_CES_ACC.pdf",width=5,height=7,useDingbats = F)
par(mar=c(5,5,1.5,0.5),mfrow=c(3,1))

plot(d_past_OSR, lwd = 2, main = "", xlab = "Operational Sex Ratio",xlim=c(0,1),
     col = "#7570B3")
mtext( "A",side = 3, adj = 0,cex=1.25)
mtext("CESM1-BGC",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(sr_past$OSR_final,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_OSR, lwd = 2,col = "#0072B2")
abline(v=mean(sr_current$OSR_final,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_ces45_OSR, lwd = 2,col = "#E6AB02")
abline(v=mean(sr_ces_45$OSR_final,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_ces85_OSR, lwd = 2,col = "#E7298A")
abline(v=mean(sr_ces_85$OSR_final,na.rm=TRUE),lty=2,col="#E7298A")
legend(0.1,6,bty="n",legend=c("Past","Present","RCP 4.5","RCP 8.5"),lwd=c(2,2,2,2),lty=1,col=c("#7570B3","#0072B2","#E6AB02","#E7298A"),cex=0.8)

plot(d_past_OSR, lwd = 2, main = "", xlab = "Operational Sex Ratio",xlim=c(0,1),
     col = "#7570B3")
mtext( "B",side = 3, adj = 0,cex=1.25)
mtext("ACCESS1-3",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(sr_past$OSR_final,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_OSR, lwd = 2,col = "#0072B2")
abline(v=mean(sr_current$OSR_final,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_acc45_OSR, lwd = 2,col = "#E6AB02")
abline(v=mean(sr_acc_45$OSR_final,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_acc85_OSR, lwd = 2,col = "#E7298A")
abline(v=mean(sr_acc_85$OSR_final,na.rm=TRUE),lty=2,col="#E7298A")


plot(d_past_OSR, lwd = 2, main = "", xlab = "Operational Sex Ratio",xlim=c(0,1),
     col = "#7570B3")
mtext( "C",side = 3, adj = 0,cex=1.25)
mtext("MIROC 5",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(sr_past$OSR_final,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_OSR, lwd = 2,col = "#0072B2")
abline(v=mean(sr_current$OSR_final,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_miroc45_OSR, lwd = 2,col = "#E6AB02")
abline(v=mean(sr_miroc_45$OSR_final,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_miroc85_OSR, lwd = 2,col = "#E7298A")
abline(v=mean(sr_miroc_85$OSR_final,na.rm=TRUE),lty=2,col="#E7298A")

dev.off()

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/POAR_SR.pdf",width=7,height=7,useDingbats = F)
par(mar=c(5,5,1.5,0.5),mfrow=c(2,2))
plot(d_past_SR, lwd = 2, main = "", xlab = "Sex ratio",xlim=c(0.48,0.52),ylim=c(0,500),
     col = "#7570B3")
mtext( "A",side = 3, adj = 0,cex=1.25)
mtext("CMCC-CM",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(sr_past$SR_final,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_SR, lwd = 2,col = "#0072B2")
abline(v=mean(sr_current$SR_final,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_cmc45_SR, lwd = 2,col = "#E6AB02")
abline(v=mean(sr_cmc_45$SR_final,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_cmc85_SR, lwd = 2,col = "#E7298A")
abline(v=mean(sr_cmc_85$SR_final,na.rm=TRUE),lty=2,col="#E7298A")


legend(0.505,450,bty="n",legend=c("Past","Present","RCP 4.5","RCP 8.5"),lwd=c(2,2,2,2),lty=1,col=c("#7570B3","#0072B2","#E6AB02","#E7298A"),cex=0.8)



plot(d_past_SR, lwd = 2, main = "", xlab = "Sex ratio",xlim=c(0.48,0.52),ylim=c(0,500),
     col = "#7570B3")
mtext( "B",side = 3, adj = 0,cex=1.25)
mtext("CESM1-BGC",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(sr_past$SR_final,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_SR, lwd = 2,col = "#0072B2")
abline(v=mean(sr_current$SR_final,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_ces45_SR, lwd = 2,col = "#E6AB02")
abline(v=mean(sr_ces_45$SR_final,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_ces85_SR, lwd = 2,col = "#E7298A")
abline(v=mean(sr_ces_85$SR_final,na.rm=TRUE),lty=2,col="#E7298A")


plot(d_past_SR, lwd = 2, main = "", xlab = "Sex ratio",xlim=c(0.48,0.52),ylim=c(0,500),
     col = "#7570B3")
mtext( "C",side = 3, adj = 0,cex=1.25)
mtext("ACCESS1-3",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(sr_past$SR_final,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_SR, lwd = 2,col = "#0072B2")
abline(v=mean(sr_current$SR_final,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_acc45_SR, lwd = 2,col = "#E6AB02")
abline(v=mean(sr_acc_45$SR_final,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_acc85_SR, lwd = 2,col = "#E7298A")
abline(v=mean(sr_acc_85$SR_final,na.rm=TRUE),lty=2,col="#E7298A")


plot(d_past_SR, lwd = 2, main = "", xlab = "Sex ratio",xlim=c(0.48,0.52),ylim=c(0,600),
     col = "#7570B3")
mtext( "D",side = 3, adj = 0,cex=1.25)
mtext("MIROC 5",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(sr_past$SR_final,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_SR, lwd = 2,col = "#0072B2")
abline(v=mean(sr_current$SR_final,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_miroc45_SR, lwd = 2,col = "#E6AB02")
abline(v=mean(sr_miroc_45$SR_final,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_miroc85_SR, lwd = 2,col = "#E7298A")
abline(v=mean(sr_miroc_85$SR_final,na.rm=TRUE),lty=2,col="#E7298A")

dev.off()

## Projection of sex ratio  across species range----
geo_sr_past <- merge(x = sr_past,y =clim_past_values,by=c("pptdorm","tempdorm"))
geo_sr_current <- merge(x = sr_current,y =clim_current_values,by=c("pptdorm","tempdorm")) 
geo_sr_cmc_85 <- merge(x = sr_cmc_85,y =climate_cmc_85_values,by=c("pptdorm","tempdorm"))
geo_sr_cmc_45 <- merge(x = sr_cmc_45,y =climate_cmc_45_values,by=c("pptdorm","tempdorm")) 
geo_sr_miroc_85 <- merge(x = sr_miroc_85,y =climate_miroc85_values,by=c("pptdorm","tempdorm"))
geo_sr_miroc_45 <- merge(x = sr_miroc_45,y =climate_miroc45_values,by=c("pptdorm","tempdorm")) 
geo_sr_ces_85 <- merge(x = sr_ces_85,y =climate_ces_85_values,by=c("pptdorm","tempdorm"))
geo_sr_ces_45 <- merge(x = sr_ces_45,y =climate_ces_45_values,by=c("pptdorm","tempdorm")) 
geo_sr_acc_85 <- merge(x = sr_acc_85 ,y =climate_acc_85_values,by=c("pptdorm","tempdorm"))
geo_sr_acc_45 <- merge(x = sr_acc_45,y =climate_acc_45_values,by=c("pptdorm","tempdorm")) 

(map_past_osr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-107,-92), ylim = c(25,41)) +
    geom_tile(data = geo_sr_past, aes(x = x, y = y, fill = OSR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    annotate(geom="text", x=-100, y=41, label="Past",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "OSR",
      colours = rainbow(10),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    labs(x = NULL, y = NULL) +
    scale_x_continuous(labels = ~ .x) +
    scale_y_continuous(labels = ~ .x) +
    theme_light()+
    theme(strip.background = element_blank(),
          legend.key.height= unit(0.5, 'cm'),
          legend.key.width= unit(0.25, 'cm'),
          legend.background=element_blank(),
          legend.title = element_text(size = 6),
          legend.text  = element_text(size = 5),
          legend.position = c(0.15, 0.72),
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="A") )

(map_current_osr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_current, aes(x = x, y = y, fill = OSR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="Current",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "OSR",
      colours = rainbow(10),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    labs(x = NULL, y = NULL) +
    scale_x_continuous(labels = ~ .x) +
    scale_y_continuous(labels = ~ .x) +
    theme_light()+
    theme(strip.background = element_blank(),
          legend.position = "none",
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7),
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="B") )

(map_cmc45_osr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_cmc_45, aes(x = x, y = y, fill = OSR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="RCP 4.5",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "OSR",
      colours = rainbow(10),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    labs(x = NULL, y = NULL) +
    scale_x_continuous(labels = ~ .x) +
    scale_y_continuous(labels = ~ .x) +
    theme_light()+
    theme(strip.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="C") )

(map_cmc85_osr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_cmc_85, aes(x = x, y = y, fill = OSR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="RCP 8.5",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "OSR",
      colours = rainbow(10),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    labs(x = NULL, y = NULL) +
    scale_x_continuous(labels = ~ .x) +
    scale_y_continuous(labels = ~ .x) +
    theme_light()+
    theme(strip.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="D") )

(map_past_sr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-107,-92), ylim = c(25,41)) +
    geom_tile(data = geo_sr_past, aes(x = x, y = y, fill = SR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    annotate(geom="text", x=-100, y=41, label="Past",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "SR",
      colours = rainbow(5),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    theme_light()+
    theme(strip.background = element_blank(),
          legend.background=element_blank(),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 5),
          legend.position = c(0.15, 0.7),
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="E") )

(map_current_sr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_current, aes(x = x, y = y, fill = SR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="Current",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "SR",
      colours = rainbow(5),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    theme_light()+
    theme(strip.background = element_blank(),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="F") )

(map_cmc45_sr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_cmc_45, aes(x = x, y = y, fill = SR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="RCP 4.5",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "SR",
      colours = rainbow(5),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    theme_light()+
    theme(strip.background = element_blank(),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="G") )
(map_cmc85_sr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_cmc_85, aes(x = x, y = y, fill = SR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="RCP 8.5",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "SR",
      colours = rainbow(5),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    theme_light()+
    theme(strip.background = element_blank(),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="H") )
Fig_geo_sr_cmc<-ggarrange(map_past_osr, map_current_osr, map_cmc45_osr, map_cmc85_osr,common.legend = FALSE,ncol = 4, nrow = 1)
ggsave("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Fig_geo_sr_cmc.pdf", Fig_geo_sr_cmc, width = 10, height = 3)



(map_ces45_osr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_ces_45, aes(x = x, y = y, fill = OSR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="RCP 4.5",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "OSR",
      colours = rainbow(5),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    labs(x = NULL, y = NULL) +
    scale_x_continuous(labels = ~ .x) +
    scale_y_continuous(labels = ~ .x) +
    theme_light()+
    theme(strip.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="C") )

(map_ces85_osr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_ces_85, aes(x = x, y = y, fill = OSR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="RCP 8.5",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "OSR",
      colours = rainbow(5),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    labs(x = NULL, y = NULL) +
    scale_x_continuous(labels = ~ .x) +
    scale_y_continuous(labels = ~ .x) +
    theme_light()+
    theme(strip.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="D") )

Fig_geo_sr_ces<-ggarrange(map_past_osr, map_current_osr, map_ces45_osr, map_ces85_osr,common.legend = FALSE,ncol = 4, nrow = 1)
ggsave("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Fig_geo_sr_ces.pdf", Fig_geo_sr_ces, width = 10, height = 4)

(map_miroc45_osr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_miroc_45, aes(x = x, y = y, fill = OSR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="RCP 4.5",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "OSR",
      colours = rainbow(5),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    labs(x = NULL, y = NULL) +
    scale_x_continuous(labels = ~ .x) +
    scale_y_continuous(labels = ~ .x) +
    theme_light()+
    theme(strip.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="C") )

(map_miroc85_osr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_miroc_85, aes(x = x, y = y, fill = OSR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="RCP 8.5",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "OSR",
      colours = rainbow(5),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    labs(x = NULL, y = NULL) +
    scale_x_continuous(labels = ~ .x) +
    scale_y_continuous(labels = ~ .x) +
    theme_light()+
    theme(strip.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="D") )

Fig_geo_sr_miroc<-ggarrange(map_past_osr, map_current_osr, map_miroc45_osr, map_miroc85_osr,common.legend = FALSE,ncol = 4, nrow = 1)
ggsave("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Fig_geo_sr_miroc.pdf", Fig_geo_sr_miroc, width = 10, height = 4)

(map_acc45_osr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_acc_45, aes(x = x, y = y, fill = OSR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="RCP 4.5",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "OSR",
      colours = rainbow(5),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    labs(x = NULL, y = NULL) +
    scale_x_continuous(labels = ~ .x) +
    scale_y_continuous(labels = ~ .x) +
    theme_light()+
    theme(strip.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="C") )

(map_acc85_osr <- ggplot()+
    # geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    # coord_sf(xlim = c(-110,-90), ylim = c(25,45)) +
    geom_tile(data = geo_sr_acc_85, aes(x = x, y = y, fill = OSR_final))+
    geom_spatvector(data =study_area,fill = NA, colour = "black")+
    # geom_point(data=poar_occ_1990_2019,aes(lon, lat),alpha=1/5)+
    annotate(geom="text", x=-100, y=41, label="RCP 8.5",
             color="grey",size = 3)+
    scale_fill_gradientn(
      name = "OSR",
      colours = rainbow(5),
      na.value = "transparent",
      breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
      limits=c(0,1))+
    labs(x = NULL, y = NULL) +
    scale_x_continuous(labels = ~ .x) +
    scale_y_continuous(labels = ~ .x) +
    theme_light()+
    theme(strip.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 7),
          legend.position = "none",
          strip.text  = element_text(face = "italic", color = "black"))+
    labs(x = "Longitude", y = "Latitude",title="D") )

Fig_geo_sr_acc<-ggarrange(map_past_osr, map_current_osr, map_acc45_osr, map_acc85_osr,common.legend = FALSE,ncol = 4, nrow = 1)
ggsave("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Fig_geo_sr_acc.pdf", Fig_geo_sr_acc, width = 10, height = 4)