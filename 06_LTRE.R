# Project: Forecasting range shifts of a dioecious plant species under climate change
# Purpose: How do sex-specific vital rates combine to determine the influence of climate variation on population growth rate (λ)? ? ( Life Table Response experiment)
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
library(ggnewscale)
library(randomForestExplainer)
library(randomForest)
library(vip)
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

data_site_pptgrow85<-data.frame(zclim_miroc85$zpptgrow,zclim_acc85$zpptgrow,zclim_cmc85$zpptgrow,zclim_ces85$zpptgrow)
data_site_tempgrow85<-data.frame(zclim_miroc85$ztempgrow,zclim_acc85$ztempgrow,zclim_cmc85$ztempgrow,zclim_ces85$ztempgrow)
data_site_pptdorm85<-data.frame(zclim_miroc85$zpptdorm,zclim_acc85$zpptdorm,zclim_cmc85$zpptdorm,zclim_ces85$zpptdorm)
data_site_tempdorm85<-data.frame(zclim_miroc85$ztempdorm,zclim_acc85$ztempdorm,zclim_cmc85$ztempdorm,zclim_ces85$ztempdorm)
data_site_85<-data.frame(zpptgrow=rowMeans(data_site_pptgrow85),ztempgrow=rowMeans(data_site_tempgrow85),zpptdorm=rowMeans(data_site_pptdorm85),ztempdorm=rowMeans(data_site_tempdorm85))



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

# Life Table Response Experiment (LTRE)-----
## Variable importance----
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
# write_csv(data_past_present_45_85,"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Forecasting Models output/Climatic data/Niche/data_past_present_45_85.csv")
# set up female and male parameter vectors to run the Matrix

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

## Estimation of population viability for two sex models and female dominant models. 
lambda_ppf<-lambda_ppf_2sex<-c()
max_yrs <- 30

for (i in 1:nrow(data_past_present_45_85)) {
  #linear model for comparison
  mat <- megamatrix_delay( F_params=F_params,
                           M_params=M_params,
                           twosex=F,
                           grow_perturb=0,
                           surv_perturb=0,
                           flow_perturb=0,
                           fert_perturb=0,
                           viab_perturb=0,
                           zpptgrow=data_past_present_45_85[i,1],
                           ztempgrow=data_past_present_45_85[i,2],
                           zpptdorm=data_past_present_45_85[i,3],
                           ztempdorm=data_past_present_45_85[i,4],
                           rfx=rfx_fun())$MEGAmat
  lambda_ppf[i] <- popbio::lambda(mat)
  # 2-sex model
  lambda_run <- lambdaSim_delay( F_params=F_params,
                                 M_params=M_params,
                                 grow_perturb=0,
                                 surv_perturb=0,
                                 flow_perturb=0,
                                 fert_perturb=0,
                                 viab_perturb=0,
                                 zpptgrow=data_past_present_45_85[i,1],
                                 ztempgrow=data_past_present_45_85[i,2],
                                 zpptdorm=data_past_present_45_85[i,3],
                                 ztempdorm=data_past_present_45_85[i,4],
                                 rfx = rfx_fun(),
                                 max.yrs = max_yrs)
  lambda_ppf_2sex[i] <- lambda_run$lambdatracker[max_yrs]
}

## Two-sex models----
Time<-c(rep("Past",14),rep("Present",14),rep("Future",28))
datltre_2sex<-data.frame(data_past_present_45_85,idorm=data_past_present_45_85$zpptdorm*data_past_present_45_85$ztempdorm,igrow=data_past_present_45_85$zpptgrow*data_past_present_45_85$ztempgrow,Time=Time,lambda=lambda_ppf_2sex)

modrf_2sex <- randomForest(lambda ~ zpptdorm + ztempdorm + zpptgrow + ztempgrow + idorm + igrow + I(zpptdorm)+I(ztempdorm)+I(zpptgrow)+I(ztempgrow), data=datltre_2sex, 
                      # three permutations per tree to estimate importance
                      importance = TRUE, nperm = 3, 
                      na.action = na.omit, mtry = 3)
plot(modrf_2sex)
v_2sex <- vip(modrf_2sex, title = "randomForest")

# ou_2sex<-tuneRF(y=datltre_2sex$lambda, x=as.data.frame(datltre_2sex[,1:6]),stepFactor = 1,improve = 0.02)
modrf_2sex<-randomForest::randomForest(y=datltre_2sex$lambda, x=as.data.frame(datltre_2sex[,1:6]),ntree=500,importance=TRUE,mtry=2)
plot(modrf_2sex$rsq,type="l")
randomForest::varImpPlot(modrf_2sex,type = 1,scale = FALSE,main="")

imp<- as.data.frame(importance(modrf))
imp$vars <- row.names(imp)
imp<-imp[,-2]
colnames(imp)<-c("Importance","Climate")
imp$Rimportance<-imp$Importance/sum(imp$Importance)
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Fig_LTRE.pdf",width=5,height=5,useDingbats = F)
(Fig_LTRE <- ggplot(data=imp, aes(x=Climate, y=Rimportance)) +
    geom_bar(stat="identity", color="black",fill="grey80", position=position_dodge())+
    labs(x="",y="Relative importance")+
    scale_x_discrete(
      breaks = c("ztempgrow","ztempdorm","zpptgrow","zpptdorm","idorm","igrow"),
      labels = c("Temperature (grow.season)", "Temperature (dorm.season)","Precipitation (grow.season)","Precipitation (dorm.season)","Precipitation*Temperature (dorm.season)","Precipitation*Temperature (grow.season)"))+
    coord_flip()+
    theme_bw()+
    theme(
      axis.text = element_text(color = "black",size = 10),
      legend.position=c(0.8, 0.2)))
dev.off()
## Female dominant model----
datltre_femal<-data.frame(data_past_present_45_85,idorm=data_past_present_45_85$zpptdorm*data_past_present_45_85$ztempdorm,igrow=data_past_present_45_85$zpptgrow*data_past_present_45_85$ztempgrow,Time=Time,lambda=lambda_ppf)


## LTRE (Male and female importance)----

pptgrow_seq<-seq(min(data_past_present_45_85$zpptgrow),max(data_past_present_45_85$zpptgrow),length.out=30)
## LTRE parameters:
# these are the indices of the intercepts and slopes
LTRE_params <- c(1,2,13,14,26,27,38,39)
# these are the effects of pptgrow on the intercepts and slopes
betas_pptgrow <- c(F_params[c("surv_pptgrow","surv_tempgrow_pptgrow","grow_pptgrow","grow_tempgrow_pptgrow","flow_pptgrow","flow_tempgrow_pptgrow","panic_pptgrow","panic_tempgrow_pptgrow")],
                   M_params[c("surv_pptgrow","surv_tempgrow_pptgrow","grow_pptgrow","grow_tempgrow_pptgrow","flow_pptgrow","flow_tempgrow_pptgrow","panic_pptgrow","panic_tempgrow_pptgrow")])
# these are the effects of pptgrow^2 on the intercepts and slopes
betas_pptgrow2 <- c(F_params[c("surv_pptgrow2","grow_pptgrow2","flow_pptgrow2","panic_pptgrow2")],
                    M_params[c("surv_pptgrow2","grow_pptgrow2","flow_pptgrow2","panic_pptgrow2")])

## Sensitivity of parameters to pptgrow -- this is just the derivative of a second-order polynomial
dp_dpptgrow_fun <- function(beta_pptgrow,beta_pptgrow2,pptgrow){return(beta_pptgrow + 2*beta_pptgrow2*pptgrow)}

## loop over pptgrow and calculated dp_dpptgrow and dlambda_dp
dp_dpptgrow <- dlambda_dp <- matrix(NA,nrow=length(LTRE_params)*2,ncol=length(pptgrow_seq))
lambda_pptgrow <- c()
perturbation <- 0.01

for(l in 1:length(pptgrow_seq)){
  lambda_pptgrow[l] <- lambdaSim_delay(F_params=F_params,M_params=M_params,zpptgrow=pptgrow_seq[l],ztempgrow=0,zpptdorm=0,ztempdorm=0,rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
  dp_dpptgrow[,l] <- dp_dpptgrow_fun(beta_pptgrow = unlist(betas_pptgrow),
                                     beta_pptgrow2 = unlist(betas_pptgrow2),
                                     pptgrow=pptgrow_seq[l])
  # cannot vectorize this part unfortunately
  for(p in 1:8){
    F_params_perturb <- F_params; M_params_perturb <- M_params;
    F_params_perturb[LTRE_params[p]] <- unlist(F_params[LTRE_params[p]]) + perturbation
    lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,zpptgrow=pptgrow_seq[l],ztempgrow=0,zpptdorm=0,ztempdorm=0,rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
    dlambda_dp[p,l] <- (lambda_perturb - lambda_pptgrow[l]) / perturbation
  }
  for(p in 9:16){
    F_params_perturb <- F_params; M_params_perturb <- M_params;
    M_params_perturb[LTRE_params[p-8]] <- unlist(M_params[LTRE_params[p-8]]) + perturbation
    lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,zpptgrow=pptgrow_seq[l],ztempgrow=0,zpptdorm=0,ztempdorm=0,rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
    dlambda_dp[p,l] <- (lambda_perturb - lambda_pptgrow[l]) / perturbation
  }
  print(l)
}
## write LTRE output based on posterior mean parameters
# write_rds(list(dp_dpptgrow=dp_dpptgrow,dlambda_dp=dlambda_dp),"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Forecasting Models output/POAR_LTREPPTGROW.rds")

## read LTRE output
POAR_LTRE_PPTGROW <- readRDS(url("https://www.dropbox.com/scl/fi/eu9xnhl0qfizaoom3uj88/POAR_LTREPPTGROW.rds?rlkey=l70j7j7n18cli9wosi4m5pgzs&dl=1"))
# put it all together 
LTRE_PPTGROW_out <- POAR_LTRE_PPTGROW$dp_dpptgrow * POAR_LTRE_PPTGROW$dlambda_dp
# Plot male and feme contribution
# ltre_cols <- c("#D55E00",NA,"#56B4E9",NA,"#009E73",NA,"#CC79A7",NA)

tempgrow_seq<-seq(min(data_past_present_45_85$ztempgrow),max(data_past_present_45_85$ztempgrow),length.out=30)
## LTRE parameters:
# these are the indices of the intercepts and slopes
LTRE_params <- c(1,2,13,14,26,27,38,39)
# these are the effects of pptgrow on the intercepts and slopes
betas_tempgrow <- c(F_params[c("surv_tempgrow","surv_tempgrow_pptgrow","grow_tempgrow","grow_tempgrow_pptgrow","flow_tempgrow","flow_tempgrow_pptgrow","panic_tempgrow","panic_tempgrow_pptgrow")],
                    M_params[c("surv_tempgrow","surv_tempgrow_pptgrow","grow_tempgrow","grow_tempgrow_pptgrow","flow_tempgrow","flow_tempgrow_pptgrow","panic_tempgrow","panic_tempgrow_pptgrow")])
# these are the effects of pptgrow^2 on the intercepts and slopes
betas_tempgrow2 <- c(F_params[c("surv_tempgrow2","grow_tempgrow2","flow_tempgrow2","panic_tempgrow2")],
                     M_params[c("surv_tempgrow2","grow_tempgrow2","flow_tempgrow2","panic_tempgrow2")])

## Sensitivity of parameters to pptgrow -- this is just the derivative of a second-order polynomial
dp_dtempgrow_fun <- function(beta_tempgrow,beta_tempgrow2,tempgrow){return(beta_tempgrow + 2*beta_tempgrow2*tempgrow)}

## loop over pptgrow and calculated dp_dpptgrow and dlambda_dp
dp_dtempgrow <- dlambda_dp <- matrix(NA,nrow=length(LTRE_params)*2,ncol=length(tempgrow_seq))
lambda_tempgrow <- c()
perturbation <- 0.01

for(l in 1:length(tempgrow_seq)){
  lambda_tempgrow[l] <- lambdaSim_delay(F_params=F_params,M_params=M_params,zpptgrow=0,ztempgrow=tempgrow_seq[l],zpptdorm=0,ztempdorm=0,rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
  dp_dtempgrow[,l] <- dp_dtempgrow_fun(beta_tempgrow = unlist(betas_tempgrow),
                                       beta_tempgrow2 = unlist(betas_tempgrow2),
                                       tempgrow=tempgrow_seq[l])
  # cannot vectorize this part unfortunately
  for(p in 1:8){
    F_params_perturb <- F_params; M_params_perturb <- M_params;
    F_params_perturb[LTRE_params[p]] <- unlist(F_params[LTRE_params[p]]) + perturbation
    lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,zpptgrow=0,ztempgrow=tempgrow_seq[l],zpptdorm=0,ztempdorm=0,rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
    dlambda_dp[p,l] <- (lambda_perturb - lambda_tempgrow[l]) / perturbation
  }
  for(p in 9:16){
    F_params_perturb <- F_params; M_params_perturb <- M_params;
    M_params_perturb[LTRE_params[p-8]] <- unlist(M_params[LTRE_params[p-8]]) + perturbation
    lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,zpptgrow=0,ztempgrow=tempgrow_seq[l],zpptdorm=0,ztempdorm=0,rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
    dlambda_dp[p,l] <- (lambda_perturb - lambda_tempgrow[l]) / perturbation
  }
  print(l)
}
## write LTRE output based on posterior mean parameters
# write_rds(list(dp_dtempgrow=dp_dtempgrow,dlambda_dp=dlambda_dp),"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Forecasting Models output/POAR_LTRETEMPGROW.rds")

## read LTRE output
POAR_LTRE_TEMPGROW <- readRDS(url("https://www.dropbox.com/scl/fi/o3zaudgef040befdqdj2j/POAR_LTRETEMPGROW.rds?rlkey=cqy74r3mw0vtvx0r6obhlmk0f&dl=1"))
# put it all together 
LTRE_TEMPGROW_out <- POAR_LTRE_TEMPGROW$dp_dtempgrow * POAR_LTRE_TEMPGROW$dlambda_dp
# Plot male and feme contribution
# ltre_cols <- c("#D55E00",NA,"#56B4E9",NA,"#009E73",NA,"#CC79A7",NA)

all_pptdorm_seq<-seq(min(data_past_present_45_85$zpptdorm),max(data_past_present_45_85$zpptdorm),length.out=30)
max_yrs<-30

## LTRE parameters:
# these are the indices of the intercepts and slopes
LTRE_params <- c(1,2,13,14,26,27,38,39)
# these are the effects of pptdorm on the intercepts and slopes
betas_pptdorm <- c(F_params[c("surv_pptdorm","surv_tempdorm_pptdorm","grow_pptdorm","grow_tempdorm_pptdorm","flow_pptdorm","flow_tempdorm_pptdorm","panic_pptdorm","panic_tempdorm_pptdorm")],
                   M_params[c("surv_pptdorm","surv_tempdorm_pptdorm","grow_pptdorm","grow_tempdorm_pptdorm","flow_pptdorm","flow_tempdorm_pptdorm","panic_pptdorm","panic_tempdorm_pptdorm")])
# these are the effects of pptdorm^2 on the intercepts and slopes
betas_pptdorm2 <- c(F_params[c("surv_pptdorm2","grow_pptdorm2","flow_pptdorm2","panic_pptdorm2")],
                    M_params[c("surv_pptdorm2","grow_pptdorm2","flow_pptdorm2","panic_pptdorm2")])

## Sensitivity of parameters to pptdorm -- this is just the derivative of a second-order polynomial
dp_dpptdorm_fun <- function(beta_pptdorm,beta_pptdorm2,pptdorm){return(beta_pptdorm + 2*beta_pptdorm2*pptdorm)}

## loop over pptdorm and calculated dp_dpptdorm and dlambda_dp
dp_dpptdorm <- dlambda_dp <- matrix(NA,nrow=length(LTRE_params)*2,ncol=length(all_pptdorm_seq))
lambda_pptdorm <- c()
perturbation <- 0.01

for(l in 1:length(all_pptdorm_seq)){
  lambda_pptdorm[l] <- lambdaSim_delay(F_params=F_params,M_params=M_params,zpptgrow=0,ztempgrow=0,ztempdorm=0,zpptdorm=all_pptdorm_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
  dp_dpptdorm[,l] <- dp_dpptdorm_fun(beta_pptdorm = unlist(betas_pptdorm),
                                     beta_pptdorm2 = unlist(betas_pptdorm2),
                                     pptdorm=all_pptdorm_seq[l])
  # cannot vectorize this part unfortunately
  for(p in 1:8){
    F_params_perturb <- F_params; M_params_perturb <- M_params;
    F_params_perturb[LTRE_params[p]] <- unlist(F_params[LTRE_params[p]]) + perturbation
    lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,zpptgrow=0,ztempgrow=0,ztempdorm=0,zpptdorm=all_pptdorm_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
    dlambda_dp[p,l] <- (lambda_perturb - lambda_pptdorm[l]) / perturbation
  }
  for(p in 9:16){
    F_params_perturb <- F_params; M_params_perturb <- M_params;
    M_params_perturb[LTRE_params[p-8]] <- unlist(M_params[LTRE_params[p-8]]) + perturbation
    lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,zpptgrow=0,ztempgrow=0,ztempdorm=0,zpptdorm=all_pptdorm_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
    dlambda_dp[p,l] <- (lambda_perturb - lambda_pptdorm[l]) / perturbation
  }
  print(l)
}
## write LTRE output based on posterior mean parameters
# write_rds(list(dp_dpptdorm=dp_dpptdorm,dlambda_dp=dlambda_dp),"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Forecasting Models output/POAR_LTREPPTDORM.rds")

POAR_LTRE_PPTDORM <- readRDS(url("https://www.dropbox.com/scl/fi/fppc6kjqx1cgayk7wf95w/POAR_LTREPPTDORM.rds?rlkey=qrv2r9atakg7e2hlhof4xoyzs&dl=1"))
# put it all together 
LTRE_PPTDORM_out <- POAR_LTRE_PPTDORM$dp_dpptdorm * POAR_LTRE_PPTDORM$dlambda_dp

all_tempdorm_seq<-seq(min(data_past_present_45_85$ztempdorm),max(data_past_present_45_85$ztempdorm),length.out=30)
max_yrs<-30

## LTRE parameters:
# these are the indices of the intercepts and slopes
LTRE_params <- c(1,2,13,14,26,27,38,39)
# these are the effects of tempdorm on the intercepts and slopes
betas_tempdorm <- c(F_params[c("surv_tempdorm","surv_tempdorm_pptdorm","grow_tempdorm","grow_tempdorm_pptdorm","flow_tempdorm","flow_tempdorm_pptdorm","panic_tempdorm","panic_tempdorm_pptdorm")],
                    M_params[c("surv_tempdorm","surv_tempdorm_pptdorm","grow_tempdorm","grow_tempdorm_pptdorm","flow_tempdorm","flow_tempdorm_pptdorm","panic_tempdorm","panic_tempdorm_pptdorm")])
# these are the effects of tempdorm^2 on the intercepts and slopes
betas_tempdorm2 <- c(F_params[c("surv_tempdorm2","grow_tempdorm2","flow_tempdorm2","panic_tempdorm2")],
                     M_params[c("surv_tempdorm2","grow_tempdorm2","flow_tempdorm2","panic_tempdorm2")])

## Sensitivity of parameters to tempdorm -- this is just the derivative of a second-order polynomial
dp_dtempdorm_fun <- function(beta_tempdorm,beta_tempdorm2,tempdorm){return(beta_tempdorm + 2*beta_tempdorm2*tempdorm)}

## loop over tempdorm and calculated dp_dtempdorm and dlambda_dp
dp_dtempdorm <- dlambda_dp <- matrix(NA,nrow=length(LTRE_params)*2,ncol=length(all_tempdorm_seq))
lambda_tempdorm <- c()
perturbation <- 0.01

for(l in 1:length(all_tempdorm_seq)){
  lambda_tempdorm[l] <- lambdaSim_delay(F_params=F_params,M_params=M_params,zpptgrow=0,ztempgrow=0,zpptdorm=0,ztempdorm=all_tempdorm_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
  dp_dtempdorm[,l] <- dp_dtempdorm_fun(beta_tempdorm = unlist(betas_tempdorm),
                                       beta_tempdorm2 = unlist(betas_tempdorm2),
                                       tempdorm=all_tempdorm_seq[l])
  # cannot vectorize this part unfortunately
  for(p in 1:8){
    F_params_perturb <- F_params; M_params_perturb <- M_params;
    F_params_perturb[LTRE_params[p]] <- unlist(F_params[LTRE_params[p]]) + perturbation
    lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,zpptgrow=0,ztempgrow=0,zpptdorm=0,ztempdorm=all_tempdorm_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
    dlambda_dp[p,l] <- (lambda_perturb - lambda_tempdorm[l]) / perturbation
  }
  for(p in 9:16){
    F_params_perturb <- F_params; M_params_perturb <- M_params;
    M_params_perturb[LTRE_params[p-8]] <- unlist(M_params[LTRE_params[p-8]]) + perturbation
    lambda_perturb <- lambdaSim_delay(F_params=F_params_perturb,M_params=M_params_perturb,zpptgrow=0,ztempgrow=0,zpptdorm=0,ztempdorm=all_tempdorm_seq[l],rfx=rfx_fun(),max.yrs=max_yrs)$lambdatracker[max_yrs]
    dlambda_dp[p,l] <- (lambda_perturb - lambda_tempdorm[l]) / perturbation
  }
  print(l)
}
## write LTRE output based on posterior mean parameters
# write_rds(list(dp_dtempdorm=dp_dtempdorm,dlambda_dp=dlambda_dp),"/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/Forecasting Models output/POAR_LTRETEMPDORM.rds")

## read LTRE output
POAR_LTRE_TEMPDORM <- readRDS(url("https://www.dropbox.com/scl/fi/l6lpwgi11lb4qkks501sh/POAR_LTRETEMPDORM.rds?rlkey=m9m0jw9qsrutr7z31rkl3hdp1&dl=1"))
# put it all together 
LTRE_TEMPDORM_out <- POAR_LTRE_TEMPDORM$dp_dtempdorm * POAR_LTRE_TEMPDORM$dlambda_dp

# ltre_cols <- c("#D55E00",NA,"#56B4E9",NA,"#009E73",NA,"#CC79A7",NA)
ltre_cols <- c("#7570B3",NA,"#0072B2",NA,"#E7298A",NA,"#FFDB6D",NA)
ltre_lty <- c(1,NA,1,NA,1,NA,1,NA)

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/LTRE_Temperature.pdf",width=8,height=7,useDingbats = F)
par(oma=c(5,3,0.5,0.5))
# par(mar=c(5,4,1,0),mfrow=c(2,2))
par(mar=c(5,5,1,1),mfrow=c(2,2))
plot(tempgrow_seq*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),colSums(LTRE_TEMPGROW_out[1:8,]),type="l",lwd=3,col="grey20",lty=1,
     xlab="Growing season temp",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"tempgrow")),ylim=c(-5,5),cex.lab=1.2)
sex.symbols(9, -3, sex = 2, cex = 1.5, lwd = 2, col = "black")  
abline(h=0,col="black")
for(i in seq(1,8,by=2)){
  lines(tempgrow_seq*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),colSums(LTRE_TEMPGROW_out[i:(i+1),]),lty=ltre_lty[i],col=ltre_cols[i],lwd=2)
}
mtext( "A",side = 3, adj = 0,cex=1.25)
# legend("topright",bty="n",legend=c("Survival","Growth","Flowering","Panicles","Total"),lwd=c(2,2,2,2,4),
#        lty=c(na.omit(ltre_lty),2),col=c(na.omit(ltre_cols),"grey20"),cex=0.8,seg.len =1)

plot(tempgrow_seq*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),colSums(LTRE_TEMPGROW_out[9:16,]),type="l",lwd=3,col="grey20",
     xlab="Growing season temp",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"tempgrow")),ylim=c(-5,5),cex.lab=1.2)
sex.symbols(9, -3, sex = 1, cex = 1.5, lwd = 2, col = "black")  
abline(h=0,col="black")
for(i in seq(9,16,by=2)){
  lines(tempgrow_seq*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),colSums(LTRE_TEMPGROW_out[i:(i+1),]),lty=ltre_lty[i-8],col=ltre_cols[i-8],lwd=2)
}
mtext( "B",side = 3, adj = 0,cex=1.25)
legend("topright",bty="n",legend=c("Survival","Growth","Flowering","Panicles","Total"),lwd=c(2,2,2,2,3),
       lty=c(na.omit(ltre_lty),1),col=c(na.omit(ltre_cols),"grey20"),cex=0.8)

# par(mar=c(5,5,1,1),mfrow=c(2,1))
plot(all_tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),colSums(LTRE_TEMPDORM_out[1:8,]),type="l",lwd=3,col="grey20",lty=1,
     xlab="Dormant season temp",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"tempdorm")),ylim=c(-5,5),cex.lab=1.2)
abline(h=0,col="black")
for(i in seq(1,8,by=2)){
  lines(all_tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),colSums(LTRE_TEMPDORM_out[i:(i+1),]),lty=ltre_lty[i],col=ltre_cols[i],lwd=2)
}
sex.symbols(25, -4, sex = 2, cex = 1.5, lwd = 2, col = "black") 
mtext( "C",side = 3, adj = 0,cex=1.25)

plot(all_tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),colSums(LTRE_TEMPDORM_out[9:16,]),type="l",lwd=3,col="grey20",
     xlab="Dormant season temp",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"tempdorm")),ylim=c(-5,5),cex.lab=1.2)
abline(h=0,col="black")
for(i in seq(9,16,by=2)){
  lines(all_tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),colSums(LTRE_TEMPDORM_out[i:(i+1),]),lty=ltre_lty[i-8],col=ltre_cols[i-8],lwd=2)
}
sex.symbols(25, -4, sex = 1, cex = 1.5, lwd = 2, col = "black") 
mtext( "D",side = 3, adj = 0,cex=1.25)
dev.off()

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/LTRE_Precipitation.pdf",width=7,height=7,useDingbats = F)
par(oma=c(5,2,1.5,0.5))
# par(mar=c(5,4,1,0),mfrow=c(2,2))
par(mar=c(5,5,1,1),mfrow=c(2,2))
plot(pptgrow_seq*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),colSums(LTRE_PPTGROW_out[1:8,]),type="l",lwd=3,ylim=c(-10,3),col="grey20",lty=1,
     xlab="Growing season precip",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"pptgrow")),cex.lab=1.2)
sex.symbols(700, 2, sex = 2, cex = 1.5, lwd = 2, col = "black")  
abline(h=0,col="black")
for(i in seq(1,8,by=2)){
  lines(pptgrow_seq*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),colSums(LTRE_PPTGROW_out[i:(i+1),]),lty=ltre_lty[i],col=ltre_cols[i],lwd=2)
}
mtext( "A",side = 3, adj = 0,cex=1.25)
plot(pptgrow_seq*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),colSums(LTRE_PPTGROW_out[9:16,]),type="l",lwd=3,col="grey20",
     xlab="Growing season precip",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"pptgrow")),cex.lab=1.2)
sex.symbols(700, 2, sex = 1, cex = 1.5, lwd = 2, col = "black")  
abline(h=0,col="black")
for(i in seq(9,16,by=2)){
  lines(pptgrow_seq*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),colSums(LTRE_PPTGROW_out[i:(i+1),]),lty=ltre_lty[i-8],col=ltre_cols[i-8],lwd=2)
}
mtext( "B",side = 3, adj = 0,cex=1.25)
legend("bottomright",bty="n",legend=c("Survival","Growth","Flowering","Panicles","Total"),lwd=c(2,2,2,2,3),
       lty=c(na.omit(ltre_lty),1),col=c(na.omit(ltre_cols),"grey20"),cex=0.8)

plot(all_pptdorm_seq*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),colSums(LTRE_PPTDORM_out[1:8,]),type="l",lwd=3,col="grey20",lty=1,
     xlab="Dormant season precip",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"pptdormant")),cex.lab=1.2)
sex.symbols(425, 2, sex = 2, cex = 1.5, lwd = 2, col = "black")  
abline(h=0,col="black")
for(i in seq(1,8,by=2)){
  lines(all_pptdorm_seq*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),colSums(LTRE_PPTDORM_out[i:(i+1),]),lty=ltre_lty[i],col=ltre_cols[i],lwd=2)
}
mtext( "C",side = 3, adj = 0,cex=1.25)

plot(all_pptdorm_seq*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),colSums(LTRE_PPTGROW_out[9:16,]),type="l",lwd=3,col="grey20",
     xlab="Dormant season precip",ylab=expression(paste(partialdiff,lambda," / ",partialdiff,"pptdorm")),cex.lab=1.2)
sex.symbols(400, -0.2, sex = 1, cex = 1.5, lwd = 2, col = "black")  
abline(h=0,col="black")
for(i in seq(9,16,by=2)){
  lines(all_pptdorm_seq*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),colSums(LTRE_PPTGROW_out[i:(i+1),]),lty=ltre_lty[i-8],col=ltre_cols[i-8],lwd=2)
}
mtext("D",side = 3, adj = 0,cex=1.25)

dev.off()