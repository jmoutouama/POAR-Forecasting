# Project: Forecasting range shifts of a dioecious plant species under climate change
# Purpose: Build Matrix Projection Models from vital (survival, growth, flowering,fertility). 
# Note: Raster files are too large to provide in public repository. They are stored on a local machine.
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
library(doParallel)
library(foreach)
library(RColorBrewer)

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
### MIROC5----
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

## Format data for past, present,future to estimate species niche across each season. 
zclim_current = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_current)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_current$zpptgrow<-(clim_current$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_current$ztempgrow<-(clim_current$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_current$zpptdorm<-(clim_current$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_current$ztempdorm<-(clim_current$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

zclim_past = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_past)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_past$zpptgrow<-(clim_past$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_past$ztempgrow<-(clim_past$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_past$zpptdorm<-(clim_past$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_past$ztempdorm<-(clim_past$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

zclim_miroc45 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_miroc45)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_miroc45$zpptgrow<-(clim_miroc45$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_miroc45$ztempgrow<-(clim_miroc45$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_miroc45$zpptdorm<-(clim_miroc45$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_miroc45$ztempdorm<-(clim_miroc45$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

zclim_acc45 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_acc45)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_acc45$zpptgrow<-(clim_access45$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_acc45$ztempgrow<-(clim_access45$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_acc45$zpptdorm<-(clim_access45$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_acc45$ztempdorm<-(clim_access45$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

zclim_cmc45 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_cmc45)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_cmc45$zpptgrow<-(clim_cmcc45$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_cmc45$ztempgrow<-(clim_cmcc45$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_cmc45$zpptdorm<-(clim_cmcc45$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_cmc45$ztempdorm<-(clim_cmcc45$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

zclim_ces45 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_ces45)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_ces45$zpptgrow<-(clim_cesm45$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_ces45$ztempgrow<-(clim_cesm45$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_ces45$zpptdorm<-(clim_cesm45$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_ces45$ztempdorm<-(clim_cesm45$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

data_site_pptgrow45<-data.frame(zclim_miroc45$zpptgrow,zclim_acc45$zpptgrow,zclim_cmc45$zpptgrow,zclim_ces45$zpptgrow)
data_site_tempgrow45<-data.frame(zclim_miroc45$ztempgrow,zclim_acc45$ztempgrow,zclim_cmc45$ztempgrow,zclim_ces45$ztempgrow)
data_site_pptdorm45<-data.frame(zclim_miroc45$zpptdorm,zclim_acc45$zpptdorm,zclim_cmc45$zpptdorm,zclim_ces45$zpptdorm)
data_site_tempdorm45<-data.frame(zclim_miroc45$ztempdorm,zclim_acc45$ztempdorm,zclim_cmc45$ztempdorm,zclim_ces45$ztempdorm)
data_site_45<-data.frame(zpptgrow=rowMeans(data_site_pptgrow45),ztempgrow=rowMeans(data_site_tempgrow45),zpptdorm=rowMeans(data_site_pptdorm45),ztempdorm=rowMeans(data_site_tempdorm45))

zclim_miroc85 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_miroc85)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_miroc85$zpptgrow<-(clim_miroc85$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_miroc85$ztempgrow<-(clim_miroc85$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_miroc85$zpptdorm<-(clim_miroc85$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_miroc85$ztempdorm<-(clim_miroc85$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

zclim_acc85 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_acc85)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_acc85$zpptgrow<-(clim_access85$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_acc85$ztempgrow<-(clim_access85$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_acc85$zpptdorm<-(clim_access85$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_acc85$ztempdorm<-(clim_access85$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

zclim_cmc85 = data.frame(matrix(nrow = 14, ncol = 4)) 
colnames(zclim_cmc85)<-c("zpptgrow","ztempgrow","zpptdorm","ztempdorm")
zclim_cmc85$zpptgrow<-(clim_cmcc85$pptgrow-mean(poar_2015_2016$pptgrow))/sd(poar_2015_2016$pptgrow)
zclim_cmc85$ztempgrow<-(clim_cmcc85$tempgrow-mean(poar_2015_2016$tempgrow))/sd(poar_2015_2016$tempgrow)
zclim_cmc85$zpptdorm<-(clim_cmcc85$pptdorm-mean(poar_2015_2016$pptdorm))/sd(poar_2015_2016$pptdorm)
zclim_cmc85$ztempdorm<-(clim_cmcc85$tempdorm-mean(poar_2015_2016$tempdorm))/sd(poar_2015_2016$tempdorm)

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
# write_csv(data_past_present_45_85,"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Forecasting Models output/Climatic data/Niche/data_past_present_45_85.csv")
data_niche<-cbind(data_past_present_45_85)

# Niche Projection----
# The estimation of the species niche (Probability of lambda being higher than 1) was estimated using a local cluster. Even with that cluster, running the code below will take a while. 
# Creating a “cluster” of CPUs 
parallel::detectCores() # found out how many cores are available to simultaneous execute different pieces of a larger computation across multiple computing processors or cores
n.cores <- parallel::detectCores() - 1 # define the number of cores we want to use
my.cluster <- parallel::makeCluster( # create the cluster
  n.cores, 
  type = "FORK" # FORK is for Mac-users. If you use a Window Machine you should use PSOCK instead. 
)
print(my.cluster) #check cluster definition (optional)
doParallel::registerDoParallel(cl = my.cluster) #register it to be used by %dopar%
foreach::getDoParRegistered() #check if it is registered (optional)
foreach::getDoParWorkers() #how many workers are available? (optional)

## Two sex models -----
### Dormant season ----
clim_dorm_post<-expand.grid(pptdorm=seq(min(data_niche$zpptdorm),max(data_niche$zpptdorm),length.out=30),tempdorm=seq(min(data_niche$ztempdorm),max(data_niche$ztempdorm),length.out=30))
n_post_draws<-300
post_draws <- sample.int(length(surv_coef$b0_s), n_post_draws)
lambda_dorm_post <-matrix(NA,nrow=n_post_draws,ncol=nrow(clim_dorm_post))
max_yrs<-30
F_params <- M_params <- list()
Timedorm<-system.time(
  lambda_dorm_post<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(clim_dorm_post),.combine='c') %dopar% {
      #set up param vectors
      ## survival
      F_params$surv_mu <- surv_coef$b0_s[post_draws[p]]
      F_params$surv_size <- surv_coef$bsize_s[post_draws[p]]
      F_params$surv_pptgrow <- surv_coef$bpptgrow_s[post_draws[p]] 
      F_params$surv_pptdorm <- surv_coef$bpptdorm_s[post_draws[p]] 
      F_params$surv_tempgrow <- surv_coef$btempgrow_s[post_draws[p]]
      F_params$surv_tempdorm <- surv_coef$btempdorm_s[post_draws[p]] 
      F_params$surv_tempgrow_pptgrow<-surv_coef$btempgrowpptgrow_s[post_draws[p]]
      F_params$surv_tempdorm_pptdorm<-surv_coef$btempdormpptdorm_s[post_draws[p]]
      F_params$surv_pptgrow2<-surv_coef$bpptgrow2_s[post_draws[p]]
      F_params$surv_pptdorm2<-surv_coef$bpptdorm2_s[post_draws[p]]
      F_params$surv_tempgrow2<-surv_coef$btempgrow2_s[post_draws[p]]
      F_params$surv_tempdorm2<-surv_coef$btempdorm2_s[post_draws[p]]
      M_params$surv_mu <- surv_coef$b0_s[post_draws[p]] + surv_coef$bsex_s[post_draws[p]]  
      M_params$surv_size <- surv_coef$bsize_s[post_draws[p]] + surv_coef$bsizesex_s[post_draws[p]]
      M_params$surv_pptgrow <- surv_coef$bpptgrow_s[post_draws[p]]  + surv_coef$bpptgrowsex_s[post_draws[p]] 
      M_params$surv_tempgrow <- surv_coef$btempgrow_s[post_draws[p]] + surv_coef$btempgrowsex_s[post_draws[p]] 
      M_params$surv_pptdorm <- surv_coef$bpptdorm_s[post_draws[p]] + surv_coef$bpptdormsex_s[post_draws[p]] 
      M_params$surv_tempdorm <- surv_coef$btempdorm_s[post_draws[p]] + surv_coef$btempdormsex_s[post_draws[p]]
      M_params$surv_tempgrow_pptgrow<-surv_coef$btempgrowpptgrow_s[post_draws[p]] + surv_coef$btempgrowpptgrowsex_s[post_draws[p]]
      M_params$surv_tempdorm_pptdorm<-surv_coef$btempdormpptdorm_s[post_draws[p]] + surv_coef$btempdormpptdormsex_s[post_draws[p]]
      M_params$surv_pptgrow2<-surv_coef$bpptgrow2_s[post_draws[p]] + surv_coef$bpptgrow2sex_s[post_draws[p]]
      M_params$surv_tempgrow2<-surv_coef$btempgrow2_s[post_draws[p]] + surv_coef$btempgrow2sex_s[post_draws[p]]
      M_params$surv_pptdorm2<-surv_coef$bpptdorm2_s[post_draws[p]] + surv_coef$bpptdorm2sex_s[post_draws[p]]
      M_params$surv_tempdorm2<-surv_coef$btempdorm2_s[post_draws[p]] + surv_coef$btempdorm2sex_s[post_draws[p]]
      ## growth
      F_params$grow_mu <- grow_coef$b0_g[post_draws[p]]
      F_params$grow_size <- grow_coef$bsize_g[post_draws[p]]
      F_params$grow_pptgrow <- grow_coef$bpptgrow_g[post_draws[p]] 
      F_params$grow_pptdorm <- grow_coef$bpptdorm_g[post_draws[p]]
      F_params$grow_tempgrow <- grow_coef$btempgrow_g[post_draws[p]] 
      F_params$grow_tempdorm <- grow_coef$btempdorm_g[post_draws[p]] 
      F_params$grow_tempgrow_pptgrow<-grow_coef$btempgrowpptgrow_g[post_draws[p]]
      F_params$grow_tempdorm_pptdorm<-grow_coef$btempdormpptdorm_g[post_draws[p]] 
      F_params$grow_pptgrow2<-grow_coef$bpptgrow2_g[post_draws[p]]
      F_params$grow_tempgrow2<-grow_coef$btempgrow2_g[post_draws[p]]
      F_params$grow_pptdorm2<-grow_coef$bpptdorm2_g[post_draws[p]] 
      F_params$grow_tempdorm2<-grow_coef$btempdorm2_g[post_draws[p]] 
      F_params$sigma_g <- grow_coef$sigma[post_draws[p]] 
      M_params$grow_mu <- grow_coef$b0_g[post_draws[p]] + grow_coef$bsex_g[post_draws[p]]  
      M_params$grow_size <- grow_coef$bsize_g[post_draws[p]] + grow_coef$bsizesex_g[post_draws[p]]
      M_params$grow_pptgrow <- grow_coef$bpptgrow_g[post_draws[p]]  + grow_coef$bpptgrowsex_g[post_draws[p]] 
      M_params$grow_pptdorm <- grow_coef$bpptdorm_g[post_draws[p]]  + grow_coef$bpptdormsex_g[post_draws[p]] 
      M_params$grow_tempgrow <- grow_coef$btempgrow_g[post_draws[p]] + grow_coef$btempgrowsex_g[post_draws[p]] 
      M_params$grow_tempdorm <- grow_coef$btempdorm_g[post_draws[p]] + grow_coef$btempdormsex_g[post_draws[p]]
      M_params$grow_tempgrow_pptgrow<-grow_coef$btempgrowpptgrow_g[post_draws[p]] + grow_coef$btempgrowpptgrowsex_g[post_draws[p]]
      M_params$grow_tempdorm_pptdorm<-grow_coef$btempdormpptdorm_g[post_draws[p]] + grow_coef$btempdormpptdormsex_g[post_draws[p]]
      M_params$grow_pptgrow2<-grow_coef$bpptgrow2_g[post_draws[p]] + grow_coef$bpptgrow2sex_g[post_draws[p]]
      M_params$grow_tempgrow2<-grow_coef$btempgrow2_g[post_draws[p]] + grow_coef$btempgrow2sex_g[post_draws[p]]
      M_params$grow_pptdorm2<-grow_coef$bpptdorm2_g[post_draws[p]] + grow_coef$bpptdorm2sex_g[post_draws[p]]
      M_params$grow_tempdorm2<-grow_coef$btempdorm2_g[post_draws[p]] + grow_coef$btempdorm2sex_g[post_draws[p]]
      M_params$sigma_g <- grow_coef$sigma[post_draws[p]] 
      ## flowering
      F_params$flow_mu <- flow_coef$b0_f[post_draws[p]]
      F_params$flow_size <- flow_coef$bsize_f[post_draws[p]]
      F_params$flow_pptgrow <- flow_coef$bpptgrow_f[post_draws[p]] 
      F_params$flow_tempgrow <- flow_coef$btempgrow_f[post_draws[p]] 
      F_params$flow_tempdorm <- flow_coef$btempdorm_f[post_draws[p]] 
      F_params$flow_pptdorm <- flow_coef$bpptdorm_f[post_draws[p]]
      F_params$flow_tempgrow_pptgrow<-flow_coef$btempgrowpptgrow_f[post_draws[p]]
      F_params$flow_tempdorm_pptdorm<-flow_coef$btempdormpptdorm_f[post_draws[p]]
      F_params$flow_pptgrow2<-flow_coef$bpptgrow2_f[post_draws[p]]
      F_params$flow_pptdorm2<-flow_coef$bpptdorm2_f[post_draws[p]]
      F_params$flow_tempgrow2<-flow_coef$btempgrow2_f[post_draws[p]]
      F_params$flow_tempdorm2<-flow_coef$btempdorm2_f[post_draws[p]]
      M_params$flow_mu <- flow_coef$b0_f[post_draws[p]] + flow_coef$bsex_f[post_draws[p]]  
      M_params$flow_size <- flow_coef$bsize_f[post_draws[p]] + flow_coef$bsizesex_f[post_draws[p]]
      M_params$flow_pptgrow <- flow_coef$bpptgrow_f[post_draws[p]]  + flow_coef$bpptgrowsex_f[post_draws[p]] 
      M_params$flow_tempgrow <- flow_coef$btempgrow_f[post_draws[p]] + flow_coef$btempgrowsex_f[post_draws[p]]
      M_params$flow_tempdorm <- flow_coef$btempdorm_f[post_draws[p]] +  flow_coef$btempdormsex_f[post_draws[p]]  
      M_params$flow_pptdorm <- flow_coef$bpptdorm_f[post_draws[p]] + flow_coef$bpptdormsex_f[post_draws[p]] 
      M_params$flow_tempgrow_pptgrow<-flow_coef$btempgrowpptgrow_f[post_draws[p]] + flow_coef$btempgrowpptgrowsex_f[post_draws[p]]
      M_params$flow_tempdorm_pptdorm<-flow_coef$btempdormpptdorm_f[post_draws[p]] + flow_coef$btempdormpptdormsex_f[post_draws[p]]
      M_params$flow_pptgrow2<-flow_coef$bpptgrow2_f[post_draws[p]] + flow_coef$bpptgrow2sex_f[post_draws[p]]
      M_params$flow_tempgrow2<-flow_coef$btempgrow2_f[post_draws[p]] + flow_coef$btempgrow2sex_f[post_draws[p]]
      M_params$flow_pptdorm2<-flow_coef$bpptdorm2_f[post_draws[p]] +  flow_coef$bpptdorm2sex_f[post_draws[p]]
      M_params$flow_tempdorm2<-flow_coef$btempdorm2_f[post_draws[p]] + flow_coef$btempdorm2sex_f[post_draws[p]]
      ## panicles
      F_params$panic_mu <- panic_coef$b0_p[post_draws[p]]
      F_params$panic_size <- panic_coef$bsize_p[post_draws[p]]
      F_params$panic_pptgrow <- panic_coef$bpptgrow_p[post_draws[p]] 
      F_params$panic_tempgrow <- panic_coef$btempgrow_p[post_draws[p]] 
      F_params$panic_tempdorm <- panic_coef$btempdorm_p[post_draws[p]] 
      F_params$panic_pptdorm <- panic_coef$bpptdorm_p[post_draws[p]]
      F_params$panic_tempgrow_pptgrow<-panic_coef$btempgrowpptgrow_p[post_draws[p]]
      F_params$panic_tempdorm_pptdorm<-panic_coef$btempdormpptdorm_p[post_draws[p]]
      F_params$panic_pptgrow2<-panic_coef$bpptgrow2_p[post_draws[p]]
      F_params$panic_tempgrow2<-panic_coef$btempgrow2_p[post_draws[p]]
      F_params$panic_pptdorm2<-panic_coef$bpptdorm2_p[post_draws[p]]
      F_params$panic_tempdorm2<-panic_coef$btempdorm2_p[post_draws[p]]
      M_params$panic_mu <- panic_coef$b0_p[post_draws[p]] + panic_coef$bsex_p[post_draws[p]]  
      M_params$panic_size <- panic_coef$bsize_p[post_draws[p]] + panic_coef$bsizesex_p[post_draws[p]]
      M_params$panic_pptgrow <- panic_coef$bpptgrow_p[post_draws[p]]  + panic_coef$bpptgrowsex_p[post_draws[p]] 
      M_params$panic_tempgrow <- panic_coef$btempgrow_p[post_draws[p]] + panic_coef$btempgrowsex_p[post_draws[p]]
      M_params$panic_tempdorm <- panic_coef$btempdorm_p[post_draws[p]] +  panic_coef$btempdormsex_p[post_draws[p]]  
      M_params$panic_pptdorm <- panic_coef$bpptdorm_p[post_draws[p]] + panic_coef$bpptdormsex_p[post_draws[p]] 
      M_params$panic_tempgrow_pptgrow<-panic_coef$btempgrowpptgrow_p[post_draws[p]] + panic_coef$btempgrowpptgrowsex_p[post_draws[p]]
      M_params$panic_tempdorm_pptdorm<-panic_coef$btempdormpptdorm_p[post_draws[p]] + panic_coef$btempdormpptdormsex_p[post_draws[p]]
      M_params$panic_pptgrow2<-panic_coef$bpptgrow2_p[post_draws[p]] + panic_coef$bpptgrow2sex_p[post_draws[p]]
      M_params$panic_tempgrow2<-panic_coef$btempgrow2_p[post_draws[p]] + panic_coef$btempgrow2sex_p[post_draws[p]]
      M_params$panic_pptdorm2<-panic_coef$bpptdorm2_p[post_draws[p]] + panic_coef$bpptdorm2sex_p[post_draws[p]]
      M_params$panic_tempdorm2<-panic_coef$btempdorm2_p[post_draws[p]] + panic_coef$btempdorm2sex_p[post_draws[p]]
      ## seed viability and misc fertility params
      F_params$v0 <- via_coef$v0[post_draws[p]] 
      F_params$a_v <- via_coef$a_v[post_draws[p]] 
      F_params$ov_per_inf <- via_coef$lambda_d[post_draws[p]] 
      F_params$germ <- via_coef$m[post_draws[p]] 
      F_params$PSR <- 0.5
      ## use POAU seedling survival for females and males
      F_params$sdlg_surv <- M_params$sdlg_surv <- sdlg_surv$sdlg_surv
      ## set max size equal between the sexes
      F_params$max_size <- M_params$max_size <- round(quantile(na.omit((poar.clim_seasonal$tillerN_t1)),probs=0.95))
      ## pull out the rfx variances
      rfx <- rfx_fun(site_tau_s = surv_coef$site_tau_s[post_draws[p]],
                     block_tau_s = surv_coef$block_tau_s[post_draws[p]],
                     source_tau_s = surv_coef$source_tau_s[post_draws[p]],
                     site_tau_g = grow_coef$site_tau_g[post_draws[p]],
                     block_tau_g = grow_coef$block_tau_g[post_draws[p]],
                     source_tau_g = grow_coef$source_tau_g[post_draws[p]],
                     site_tau_f = flow_coef$site_tau_f[post_draws[p]],
                     block_tau_f = flow_coef$block_tau_f[post_draws[p]],
                     source_tau_f = flow_coef$source_tau_f[post_draws[p]],
                     site_tau_p = panic_coef$site_tau_p[post_draws[p]],
                     block_tau_p = panic_coef$block_tau_p[post_draws[p]],
                     source_tau_p = panic_coef$source_tau_p[post_draws[p]])
      #Estimation+process error
      lambdaSim_delay(F_params=F_params,
                      M_params=M_params,
                      grow_perturb=0,
                      surv_perturb=0,
                      flow_perturb=0,
                      fert_perturb=0,
                      viab_perturb=0,
                      zpptgrow=0,
                      ztempgrow=0,
                      zpptdorm=clim_dorm_post[l,1],
                      ztempdorm=clim_dorm_post[l,2],
                      rfx = rfx_fun(),
                      max.yrs = max_yrs)$lambdatracker[max_yrs]
      
    }
)

## write output as a dataframe
data_dormant_niche<-data.frame(lambda_dorm_post)
# write_rds(data_dormant_niche,"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Forecasting Models output/data_dormant_niche.rds")
nichedormant<-readRDS(url("https://www.dropbox.com/scl/fi/ejueiyu58g96eqlb60js4/data_dormant_niche.rds?rlkey=8xjqn6reelmqtcbzm1na61g2m&dl=1"))
prob_lambda_dorm_post<-c()
for(l in 1:ncol(nichedormant)){
  prob_lambda_dorm_post[l]<-mean(nichedormant[,l]>1,na.rm=T)
}
all_pptdorm_seq<-seq(min(data_niche$zpptdorm*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm)),max(data_niche$zpptdorm*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm)),length.out=30)
all_tempdorm_seq<-seq(min(data_niche$ztempdorm*sd(poar_2015_2016$tempdorm)+mean(poar_2015_2016$tempdorm)),max(data_niche$ztempdorm*sd(poar_2015_2016$tempdorm)+mean(poar_2015_2016$tempdorm)),length.out=30)

mtrx_dorm_2sex <- matrix(prob_lambda_dorm_post, nrow = 30, dimnames = list(all_pptdorm_seq, all_tempdorm_seq))

### Growing season ----
clim_grow_post<-expand.grid(pptgrow=seq(min(data_niche$zpptgrow),max(data_niche$zpptgrow),length.out=30),tempgrow=seq(min(data_niche$ztempgrow),max(data_niche$ztempgrow),length.out=30))
n_post_draws<-300
post_draws <- sample.int(length(surv_coef$b0_s), n_post_draws)
lambda_grow_post <-matrix(NA,nrow=n_post_draws,ncol=nrow(clim_grow_post))
max_yrs<-30
F_params <- M_params <- list()
Timegrow<-system.time(
  lambda_grow_post<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(clim_grow_post),.combine='c') %dopar% {
      #set up param vectors
      ## survival
      F_params$surv_mu <- surv_coef$b0_s[post_draws[p]]
      F_params$surv_size <- surv_coef$bsize_s[post_draws[p]]
      F_params$surv_pptgrow <- surv_coef$bpptgrow_s[post_draws[p]] 
      F_params$surv_pptdorm <- surv_coef$bpptdorm_s[post_draws[p]] 
      F_params$surv_tempgrow <- surv_coef$btempgrow_s[post_draws[p]]
      F_params$surv_tempdorm <- surv_coef$btempdorm_s[post_draws[p]] 
      F_params$surv_tempgrow_pptgrow<-surv_coef$btempgrowpptgrow_s[post_draws[p]]
      F_params$surv_tempdorm_pptdorm<-surv_coef$btempdormpptdorm_s[post_draws[p]]
      F_params$surv_pptgrow2<-surv_coef$bpptgrow2_s[post_draws[p]]
      F_params$surv_pptdorm2<-surv_coef$bpptdorm2_s[post_draws[p]]
      F_params$surv_tempgrow2<-surv_coef$btempgrow2_s[post_draws[p]]
      F_params$surv_tempdorm2<-surv_coef$btempdorm2_s[post_draws[p]]
      M_params$surv_mu <- surv_coef$b0_s[post_draws[p]] + surv_coef$bsex_s[post_draws[p]]  
      M_params$surv_size <- surv_coef$bsize_s[post_draws[p]] + surv_coef$bsizesex_s[post_draws[p]]
      M_params$surv_pptgrow <- surv_coef$bpptgrow_s[post_draws[p]]  + surv_coef$bpptgrowsex_s[post_draws[p]] 
      M_params$surv_tempgrow <- surv_coef$btempgrow_s[post_draws[p]] + surv_coef$btempgrowsex_s[post_draws[p]] 
      M_params$surv_pptdorm <- surv_coef$bpptdorm_s[post_draws[p]] + surv_coef$bpptdormsex_s[post_draws[p]] 
      M_params$surv_tempdorm <- surv_coef$btempdorm_s[post_draws[p]] + surv_coef$btempdormsex_s[post_draws[p]]
      M_params$surv_tempgrow_pptgrow<-surv_coef$btempgrowpptgrow_s[post_draws[p]] + surv_coef$btempgrowpptgrowsex_s[post_draws[p]]
      M_params$surv_tempdorm_pptdorm<-surv_coef$btempdormpptdorm_s[post_draws[p]] + surv_coef$btempdormpptdormsex_s[post_draws[p]]
      M_params$surv_pptgrow2<-surv_coef$bpptgrow2_s[post_draws[p]] + surv_coef$bpptgrow2sex_s[post_draws[p]]
      M_params$surv_tempgrow2<-surv_coef$btempgrow2_s[post_draws[p]] + surv_coef$btempgrow2sex_s[post_draws[p]]
      M_params$surv_pptdorm2<-surv_coef$bpptdorm2_s[post_draws[p]] + surv_coef$bpptdorm2sex_s[post_draws[p]]
      M_params$surv_tempdorm2<-surv_coef$btempdorm2_s[post_draws[p]] + surv_coef$btempdorm2sex_s[post_draws[p]]
      ## growth
      F_params$grow_mu <- grow_coef$b0_g[post_draws[p]]
      F_params$grow_size <- grow_coef$bsize_g[post_draws[p]]
      F_params$grow_pptgrow <- grow_coef$bpptgrow_g[post_draws[p]] 
      F_params$grow_pptdorm <- grow_coef$bpptdorm_g[post_draws[p]]
      F_params$grow_tempgrow <- grow_coef$btempgrow_g[post_draws[p]] 
      F_params$grow_tempdorm <- grow_coef$btempdorm_g[post_draws[p]] 
      F_params$grow_tempgrow_pptgrow<-grow_coef$btempgrowpptgrow_g[post_draws[p]]
      F_params$grow_tempdorm_pptdorm<-grow_coef$btempdormpptdorm_g[post_draws[p]] 
      F_params$grow_pptgrow2<-grow_coef$bpptgrow2_g[post_draws[p]]
      F_params$grow_tempgrow2<-grow_coef$btempgrow2_g[post_draws[p]]
      F_params$grow_pptdorm2<-grow_coef$bpptdorm2_g[post_draws[p]] 
      F_params$grow_tempdorm2<-grow_coef$btempdorm2_g[post_draws[p]] 
      F_params$sigma_g <- grow_coef$sigma[post_draws[p]] 
      M_params$grow_mu <- grow_coef$b0_g[post_draws[p]] + grow_coef$bsex_g[post_draws[p]]  
      M_params$grow_size <- grow_coef$bsize_g[post_draws[p]] + grow_coef$bsizesex_g[post_draws[p]]
      M_params$grow_pptgrow <- grow_coef$bpptgrow_g[post_draws[p]]  + grow_coef$bpptgrowsex_g[post_draws[p]] 
      M_params$grow_pptdorm <- grow_coef$bpptdorm_g[post_draws[p]]  + grow_coef$bpptdormsex_g[post_draws[p]] 
      M_params$grow_tempgrow <- grow_coef$btempgrow_g[post_draws[p]] + grow_coef$btempgrowsex_g[post_draws[p]] 
      M_params$grow_tempdorm <- grow_coef$btempdorm_g[post_draws[p]] + grow_coef$btempdormsex_g[post_draws[p]]
      M_params$grow_tempgrow_pptgrow<-grow_coef$btempgrowpptgrow_g[post_draws[p]] + grow_coef$btempgrowpptgrowsex_g[post_draws[p]]
      M_params$grow_tempdorm_pptdorm<-grow_coef$btempdormpptdorm_g[post_draws[p]] + grow_coef$btempdormpptdormsex_g[post_draws[p]]
      M_params$grow_pptgrow2<-grow_coef$bpptgrow2_g[post_draws[p]] + grow_coef$bpptgrow2sex_g[post_draws[p]]
      M_params$grow_tempgrow2<-grow_coef$btempgrow2_g[post_draws[p]] + grow_coef$btempgrow2sex_g[post_draws[p]]
      M_params$grow_pptdorm2<-grow_coef$bpptdorm2_g[post_draws[p]] + grow_coef$bpptdorm2sex_g[post_draws[p]]
      M_params$grow_tempdorm2<-grow_coef$btempdorm2_g[post_draws[p]] + grow_coef$btempdorm2sex_g[post_draws[p]]
      M_params$sigma_g <- grow_coef$sigma[post_draws[p]] 
      ## flowering
      F_params$flow_mu <- flow_coef$b0_f[post_draws[p]]
      F_params$flow_size <- flow_coef$bsize_f[post_draws[p]]
      F_params$flow_pptgrow <- flow_coef$bpptgrow_f[post_draws[p]] 
      F_params$flow_tempgrow <- flow_coef$btempgrow_f[post_draws[p]] 
      F_params$flow_tempdorm <- flow_coef$btempdorm_f[post_draws[p]] 
      F_params$flow_pptdorm <- flow_coef$bpptdorm_f[post_draws[p]]
      F_params$flow_tempgrow_pptgrow<-flow_coef$btempgrowpptgrow_f[post_draws[p]]
      F_params$flow_tempdorm_pptdorm<-flow_coef$btempdormpptdorm_f[post_draws[p]]
      F_params$flow_pptgrow2<-flow_coef$bpptgrow2_f[post_draws[p]]
      F_params$flow_pptdorm2<-flow_coef$bpptdorm2_f[post_draws[p]]
      F_params$flow_tempgrow2<-flow_coef$btempgrow2_f[post_draws[p]]
      F_params$flow_tempdorm2<-flow_coef$btempdorm2_f[post_draws[p]]
      M_params$flow_mu <- flow_coef$b0_f[post_draws[p]] + flow_coef$bsex_f[post_draws[p]]  
      M_params$flow_size <- flow_coef$bsize_f[post_draws[p]] + flow_coef$bsizesex_f[post_draws[p]]
      M_params$flow_pptgrow <- flow_coef$bpptgrow_f[post_draws[p]]  + flow_coef$bpptgrowsex_f[post_draws[p]] 
      M_params$flow_tempgrow <- flow_coef$btempgrow_f[post_draws[p]] + flow_coef$btempgrowsex_f[post_draws[p]]
      M_params$flow_tempdorm <- flow_coef$btempdorm_f[post_draws[p]] +  flow_coef$btempdormsex_f[post_draws[p]]  
      M_params$flow_pptdorm <- flow_coef$bpptdorm_f[post_draws[p]] + flow_coef$bpptdormsex_f[post_draws[p]] 
      M_params$flow_tempgrow_pptgrow<-flow_coef$btempgrowpptgrow_f[post_draws[p]] + flow_coef$btempgrowpptgrowsex_f[post_draws[p]]
      M_params$flow_tempdorm_pptdorm<-flow_coef$btempdormpptdorm_f[post_draws[p]] + flow_coef$btempdormpptdormsex_f[post_draws[p]]
      M_params$flow_pptgrow2<-flow_coef$bpptgrow2_f[post_draws[p]] + flow_coef$bpptgrow2sex_f[post_draws[p]]
      M_params$flow_tempgrow2<-flow_coef$btempgrow2_f[post_draws[p]] + flow_coef$btempgrow2sex_f[post_draws[p]]
      M_params$flow_pptdorm2<-flow_coef$bpptdorm2_f[post_draws[p]] +  flow_coef$bpptdorm2sex_f[post_draws[p]]
      M_params$flow_tempdorm2<-flow_coef$btempdorm2_f[post_draws[p]] + flow_coef$btempdorm2sex_f[post_draws[p]]
      ## panicles
      F_params$panic_mu <- panic_coef$b0_p[post_draws[p]]
      F_params$panic_size <- panic_coef$bsize_p[post_draws[p]]
      F_params$panic_pptgrow <- panic_coef$bpptgrow_p[post_draws[p]] 
      F_params$panic_tempgrow <- panic_coef$btempgrow_p[post_draws[p]] 
      F_params$panic_tempdorm <- panic_coef$btempdorm_p[post_draws[p]] 
      F_params$panic_pptdorm <- panic_coef$bpptdorm_p[post_draws[p]]
      F_params$panic_tempgrow_pptgrow<-panic_coef$btempgrowpptgrow_p[post_draws[p]]
      F_params$panic_tempdorm_pptdorm<-panic_coef$btempdormpptdorm_p[post_draws[p]]
      F_params$panic_pptgrow2<-panic_coef$bpptgrow2_p[post_draws[p]]
      F_params$panic_tempgrow2<-panic_coef$btempgrow2_p[post_draws[p]]
      F_params$panic_pptdorm2<-panic_coef$bpptdorm2_p[post_draws[p]]
      F_params$panic_tempdorm2<-panic_coef$btempdorm2_p[post_draws[p]]
      M_params$panic_mu <- panic_coef$b0_p[post_draws[p]] + panic_coef$bsex_p[post_draws[p]]  
      M_params$panic_size <- panic_coef$bsize_p[post_draws[p]] + panic_coef$bsizesex_p[post_draws[p]]
      M_params$panic_pptgrow <- panic_coef$bpptgrow_p[post_draws[p]]  + panic_coef$bpptgrowsex_p[post_draws[p]] 
      M_params$panic_tempgrow <- panic_coef$btempgrow_p[post_draws[p]] + panic_coef$btempgrowsex_p[post_draws[p]]
      M_params$panic_tempdorm <- panic_coef$btempdorm_p[post_draws[p]] +  panic_coef$btempdormsex_p[post_draws[p]]  
      M_params$panic_pptdorm <- panic_coef$bpptdorm_p[post_draws[p]] + panic_coef$bpptdormsex_p[post_draws[p]] 
      M_params$panic_tempgrow_pptgrow<-panic_coef$btempgrowpptgrow_p[post_draws[p]] + panic_coef$btempgrowpptgrowsex_p[post_draws[p]]
      M_params$panic_tempdorm_pptdorm<-panic_coef$btempdormpptdorm_p[post_draws[p]] + panic_coef$btempdormpptdormsex_p[post_draws[p]]
      M_params$panic_pptgrow2<-panic_coef$bpptgrow2_p[post_draws[p]] + panic_coef$bpptgrow2sex_p[post_draws[p]]
      M_params$panic_tempgrow2<-panic_coef$btempgrow2_p[post_draws[p]] + panic_coef$btempgrow2sex_p[post_draws[p]]
      M_params$panic_pptdorm2<-panic_coef$bpptdorm2_p[post_draws[p]] + panic_coef$bpptdorm2sex_p[post_draws[p]]
      M_params$panic_tempdorm2<-panic_coef$btempdorm2_p[post_draws[p]] + panic_coef$btempdorm2sex_p[post_draws[p]]
      ## seed viability and misc fertility params
      F_params$v0 <- via_coef$v0[post_draws[p]] 
      F_params$a_v <- via_coef$a_v[post_draws[p]] 
      F_params$ov_per_inf <- via_coef$lambda_d[post_draws[p]] 
      F_params$germ <- via_coef$m[post_draws[p]] 
      F_params$PSR <- 0.5
      ## use POAU seedling survival for females and males
      F_params$sdlg_surv <- M_params$sdlg_surv <- sdlg_surv$sdlg_surv
      ## set max size equal between the sexes
      F_params$max_size <- M_params$max_size <- round(quantile(na.omit((poar.clim_seasonal$tillerN_t1)),probs=0.95))
      ## pull out the rfx variances
      rfx <- rfx_fun(site_tau_s = surv_coef$site_tau_s[post_draws[p]],
                     block_tau_s = surv_coef$block_tau_s[post_draws[p]],
                     source_tau_s = surv_coef$source_tau_s[post_draws[p]],
                     site_tau_g = grow_coef$site_tau_g[post_draws[p]],
                     block_tau_g = grow_coef$block_tau_g[post_draws[p]],
                     source_tau_g = grow_coef$source_tau_g[post_draws[p]],
                     site_tau_f = flow_coef$site_tau_f[post_draws[p]],
                     block_tau_f = flow_coef$block_tau_f[post_draws[p]],
                     source_tau_f = flow_coef$source_tau_f[post_draws[p]],
                     site_tau_p = panic_coef$site_tau_p[post_draws[p]],
                     block_tau_p = panic_coef$block_tau_p[post_draws[p]],
                     source_tau_p = panic_coef$source_tau_p[post_draws[p]])
      lambdaSim_delay(F_params=F_params,
                      M_params=M_params,
                      grow_perturb=0,
                      surv_perturb=0,
                      flow_perturb=0,
                      fert_perturb=0,
                      viab_perturb=0,
                      zpptgrow=clim_grow_post[l,1],
                      ztempgrow=clim_grow_post[l,2],
                      zpptdorm=0,
                      ztempdorm=0,
                      rfx = rfx_fun(),
                      max.yrs = max_yrs)$lambdatracker[max_yrs]
      # mean(lambda_grow_post[,l]>1,na.rm=T)
      
    }
)

data_growing_niche<-data.frame(lambda_grow_post)
## write output as a rds and csv and rds
# write_rds(data_growing_niche,"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Forecasting Models output/data_growing_niche.rds")
nichegrowing <- readRDS(url("https://www.dropbox.com/scl/fi/ggxrqzgz2wwmo1mec2d5m/data_growing_niche.rds?rlkey=6uv4l0e2ajovdit7gt778aexw&dl=1"))
prob_lambda_grow_post<-c()
for(l in 1:ncol(nichegrowing)){
  prob_lambda_grow_post[l]<-mean(nichegrowing[,l]>1,na.rm=T)
}
all_pptgrow_seq<-seq(min(data_niche$zpptgrow*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow)),max(data_niche$zpptgrow*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow)),length.out=30)
all_tempgrow_seq<-seq(min(data_niche$ztempgrow*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow)),max(data_niche$ztempgrow*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow)),length.out=30)
mtrx_grow_2sex <- matrix(prob_lambda_grow_post, nrow = 30, dimnames = list(all_pptgrow_seq,all_tempgrow_seq))

# Setting climate values 9past, present and future) for each site. 

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

## Female dominant models ----
### Dormant season -----
nichedormant_fd<-readRDS(url("https://www.dropbox.com/scl/fi/u0yeohxrienfp96hopvl5/data_dormant_niche.rds?rlkey=lhum5y9ymqtpkegkw75seruuc&dl=1"))
prob_lambda_dorm_post_fd<-c()
for(l in 1:ncol(nichedormant_fd)){
  prob_lambda_dorm_post_fd[l]<-mean(nichedormant_fd[,l]>1,na.rm=T)
}

all_pptdorm_seq<-seq(min(data_niche$zpptdorm*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm)),max(data_niche$zpptdorm*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm)),length.out=30)
all_tempdorm_seq<-seq(min(data_niche$ztempdorm*sd(poar_2015_2016$tempdorm)+mean(poar_2015_2016$tempdorm)),max(data_niche$ztempdorm*sd(poar_2015_2016$tempdorm)+mean(poar_2015_2016$tempdorm)),length.out=30)

mtrx_dorm_fd <- matrix(prob_lambda_dorm_post_fd, nrow = 30, dimnames = list(all_pptdorm_seq,all_tempdorm_seq))

### Growing season -----
nichegrowing_fd <- readRDS(url("https://www.dropbox.com/scl/fi/xgi4mzlgwrw1k73n6sdkw/data_growing_niche.rds?rlkey=7gw2nh7mw67ind2dhs9mbfszy&dl=1"))
prob_lambda_grow_post_fd<-c()
for(l in 1:ncol(nichegrowing_fd)){
  prob_lambda_grow_post_fd[l]<-mean(nichegrowing_fd[,l]>1,na.rm=T)
}
all_pptgrow_seq<-seq(min(data_niche$zpptgrow*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow)),max(data_niche$zpptgrow*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow)),length.out=30)
all_tempgrow_seq<-seq(min(data_niche$ztempgrow*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow)),max(data_niche$ztempgrow*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow)),length.out=30)

mtrx_grow_fd <- matrix(prob_lambda_grow_post_fd, nrow = 30, dimnames = list(all_pptgrow_seq,all_tempgrow_seq))

## Difference between Female dominant and two sex models

niche_diff_dormant<-prob_lambda_dorm_post_fd-prob_lambda_dorm_post
mtrx_diff_dorm <- matrix(niche_diff_dormant, nrow = 30, dimnames = list(all_pptdorm_seq,all_tempdorm_seq))
niche_diff_grow<-prob_lambda_grow_post_fd-prob_lambda_grow_post
mtrx_diff_grow<- matrix(niche_diff_grow, nrow = 30, dimnames = list(all_pptgrow_seq,all_tempgrow_seq))

# Fig 3. Niche 

scal_breaks <- c(seq(0, 1, length.out = 101))
# scal_breaks <- seq(0, 1, length.out=101,labels=c(0.00, 0.25, 0.50, 0.75, 1.00),limits=c(0,1)
scal_breaks_diff <- c(seq(-0.2, 0.80, length.out = 101))
scal_breaks_diff_dorm <- c(seq(-0.05, 0.30, length.out = 101))

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/niche_dormant_growing.pdf",width=8,height=10,useDingbats = F)
par(mar=c(5,5,2,2),mfrow=c(3,2))
fields::image.plot(all_pptdorm_seq,all_tempdorm_seq,mtrx_dorm_fd,col=topo.colors(100),xlab="Dormant season precip",ylab="Dormant season temp",main="",breaks=scal_breaks,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste("Pr (",lambda, ">1)")), side=3, font=2, line=0.3, cex=0.5))  
# contour(all_pptdorm_seq,all_tempdorm_seq,mtrx_dorm_fd,add=T,labcex=0.75,col="black") 
mtext( "A",side = 3, adj = 0,cex=1.2)
mtext("Female-dominant model (F)",side = 3, adj = 0.5,cex=1,line=0.3)
text(x=x_current_d,y=y_current_d,label = paste0("+"),col="black",cex=2)
text(x=x_past_d,y=y_past_d,label = paste0("o"),col="black",cex=2)
text(x=x_miroc45_d,y=y_miroc45_d,label = paste0("*"),col="black",cex=2)
text(x=x_miroc85_d,y=y_miroc85_d,label = paste0("-"),col="black",cex=2)
fields::image.plot(all_pptgrow_seq,all_tempgrow_seq,mtrx_grow_fd,col=topo.colors(100),xlab="Growing season precip",ylab="Growing season tem",main="",breaks=scal_breaks,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste("Pr (",lambda, ">1)")), side=3, font=2, line=0.3, cex=0.5))  
# contour(all_pptgrow_seq,all_tempgrow_seq,mtrx_grow_fd,add=T,labcex=0.75,col="black") 
mtext( "B",side = 3, adj = 0,cex=1.2)
mtext("Female-dominant model (F)",side = 3, adj = 0.5,cex=1,line=0.3)
text(x=x_current_g,y=y_current_g,label = paste0("+"),col="black",cex=2)
text(x=x_past_g,y=y_past_g,label = paste0("o"),col="black",cex=2)
text(x=x_miroc45_g,y=y_miroc45_g,label = paste0("*"),col="black",cex=2)
text(x=x_miroc85_g,y=y_miroc85_g,label = paste0("-"),col="black",cex=2)
fields::image.plot(all_pptdorm_seq,all_tempdorm_seq,mtrx_dorm_2sex,col=topo.colors(100),xlab="Dormant season precip",ylab="Dormant season temp",main="",breaks=scal_breaks,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste("Pr (",lambda, ">1)")), side=3, font=2, line=0.3, cex=0.5))  
# contour(all_pptdorm_seq,all_tempdorm_seq,mtrx_dorm_2sex,add=T,labcex=0.75,col="black") 
mtext( "C",side = 3, adj = 0,cex=1.2)
mtext("Two-sex model (FM)",side = 3, adj = 0.5,cex=1,line=0.3)
text(x=x_current_d,y=y_current_d,label = paste0("+"),col="black",cex=2)
text(x=x_past_d,y=y_past_d,label = paste0("o"),col="black",cex=2)
text(x=x_miroc45_d,y=y_miroc45_d,label = paste0("*"),col="black",cex=2)
text(x=x_miroc85_d,y=y_miroc85_d,label = paste0("-"),col="black",cex=2)
fields::image.plot(all_pptgrow_seq,all_tempgrow_seq,mtrx_grow_2sex,col=topo.colors(100),xlab="Growing season precip",ylab="Growing season temp",main="",breaks=scal_breaks,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste("Pr (",lambda, ">1)")), side=3, font=2, line=0.3, cex=0.5))  
# contour(all_pptgrow_seq,all_tempgrow_seq,mtrx_grow_2sex,add=T,labcex=0.75,col="black") 
mtext( "D",side = 3, adj = 0,cex=1.2)
mtext("Two-sex model (FM)",side = 3, adj = 0.5,cex=1,line=0.3)
text(x=x_current_g,y=y_current_g,label = paste0("+"),col="black",cex=2)
text(x=x_past_g,y=y_past_g,label = paste0("o"),col="black",cex=2)
text(x=x_miroc45_g,y=y_miroc45_g,label = paste0("*"),col="black",cex=2)
text(x=x_miroc85_g,y=y_miroc85_g,label = paste0("-"),col="black",cex=2)
fields::image.plot(all_pptdorm_seq,all_tempdorm_seq,mtrx_diff_dorm,col=colorspace::diverge_hcl(100),xlab="Dormant season precip",ylab="Dormant season temp",main="" ,breaks=scal_breaks_diff_dorm,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste("Pr (",lambda, ">1)")), side=3, font=2, line=0.3, cex=0.5))  
# contour(all_pptdorm_seq,all_tempdorm_seq,mtrx_diff_dorm,add=T,labcex=0.75,col="black") 
mtext( "E",side = 3, adj = 0,cex=1.2)
mtext("F-FM",side = 3, adj = 0.5,cex=1,line=0.3)
text(x=x_current_d,y=y_current_d,label = paste0("+"),col="black",cex=2)
text(x=x_past_d,y=y_past_d,label = paste0("o"),col="black",cex=2)
text(x=x_miroc45_d,y=y_miroc45_d,label = paste0("*"),col="black",cex=2)
text(x=x_miroc85_d,y=y_miroc85_d,label = paste0("-"),col="black",cex=2)
fields::image.plot(all_pptgrow_seq,all_tempgrow_seq,mtrx_diff_grow,col=colorspace::diverge_hcl(100),xlab="Growing season precip",ylab="Growing season temp",main="",breaks=scal_breaks_diff,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste("Pr (",lambda, ">1)")), side=3, font=2, line=0.3, cex=0.5))  
# contour(all_pptgrow_seq,all_tempgrow_seq,mtrx_diff_grow,add=T,labcex=0.75,col="black") 
mtext( "F",side = 3, adj = 0,cex=1.2)
mtext("F-FM",side = 3, adj = 0.5,cex=1,line=0.3)
text(x=x_current_g,y=y_current_g,label = paste0("+"),col="black",cex=2)
text(x=x_past_g,y=y_past_g,label = paste0("o"),col="black",cex=2)
text(x=x_miroc45_g,y=y_miroc45_g,label = paste0("*"),col="black",cex=2)
text(x=x_miroc85_g,y=y_miroc85_g,label = paste0("-"),col="black",cex=2)
dev.off()


