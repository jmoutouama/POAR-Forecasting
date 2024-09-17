# Project: Forecasting range shifts of a dioecious plant species under climate change
# Purpose: Build Matrix Projection Models from vital (survival, growth, flowering,fertility). 
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
Conditions_cmcc_85 <- terra::extract(clim_cmcc_85,Poa_longlat)
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

# Relation lambda and each climatic variable -----
cbPalette <- c("#999999", "#E69F00")
```
