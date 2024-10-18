# Project: Forecasting range shifts of a dioecious plant species under climate change
# Question: How do sex-specific vital rates combine to determine the influence of climate variation on population growth rate (λ) ?
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

# Vital rates information ----
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

## Format data for past, present,future to project population viability across space
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

# Population viability  Projection----
# The estimation of the species niche (Probability of lambda being higher than 1) was estimated using a local cluster. Even with that cluster, running the code below will take a while (weeks depending on number of core available on you local machine).
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

## Female dominant models ----
n_post_draws<-100 
post_draws <- sample.int(length(surv_coef$b0_s), n_post_draws)

### Present----
# lambda_post_dorm <- Fdom_lambda_dorm  <- matrix(NA,nrow=n_post_draws,ncol=length(pptdorm_seq))
lambda_post_current<-matrix(NA,nrow=n_post_draws,ncol=nrow(clim_current_values))
dim(lambda_post_current)
# lambda_dorm_post<-c()

F_params <- M_params <- list()
Timecurrent<-system.time(
  lambda_post_current<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(clim_current_values),.combine='c') %dopar% {
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
      mat_current_post_rfx <- megamatrix_delay(F_params=F_params,
                                               M_params=M_params,
                                               twosex=F,
                                               grow_perturb=0,
                                               surv_perturb=0,
                                               flow_perturb=0,
                                               fert_perturb=0,
                                               viab_perturb=0,
                                               zpptgrow=clim_current_values[l,5],
                                               ztempgrow=clim_current_values[l,6],
                                               zpptdorm=clim_current_values[l,7],
                                               ztempdorm=clim_current_values[l,8],
                                               rfx=rfx_fun())$MEGAmat
      popbio::lambda(mat_current_post_rfx)
      
    }
)

geo_lambbacurrent_fd<-data.frame(lambda_post_current)

### Past----
clim_past_values_clean<-na.omit(clim_past_values)
summary(clim_past_values_clean)
dim(clim_past_values_clean)
n_post_draws<-100 
post_draws <- sample.int(length(surv_coef$b0_s), n_post_draws)
# lambda_post_dorm <- Fdom_lambda_dorm  <- matrix(NA,nrow=n_post_draws,ncol=length(pptdorm_seq))
lambda_post_past<-matrix(NA,nrow=n_post_draws,ncol=nrow(clim_past_values_clean))
dim(lambda_post_past)
# lambda_dorm_post<-c()

F_params <- M_params <- list()
Timepast<-system.time(
  lambda_post_past<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(clim_past_values_clean),.combine='c') %dopar% {
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
      mat_past_post_rfx <- megamatrix_delay(F_params=F_params,
                                            M_params=M_params,
                                            twosex=F,
                                            grow_perturb=0,
                                            surv_perturb=0,
                                            flow_perturb=0,
                                            fert_perturb=0,
                                            viab_perturb=0,
                                            zpptgrow=clim_past_values_clean[l,5],
                                            ztempgrow=clim_past_values_clean[l,6],
                                            zpptdorm=clim_past_values_clean[l,7],
                                            ztempdorm=clim_past_values_clean[l,8],
                                            rfx=rfx_fun())$MEGAmat
      popbio::lambda(mat_past_post_rfx)
      
    }
)

geo_lambbapast_fd <-data.frame(lambda_post_past)

### Future----
#### ACCESS 45----
lambda_post_acc45<-matrix(NA,nrow=n_post_draws,ncol=nrow(climate_acc_45_values))
dim(lambda_post_acc45)
# lambda_dorm_post<-c()

F_params <- M_params <- list()
Timeacc45<-system.time(
  lambda_acc45<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(climate_acc_45_values),.combine='c') %dopar% {
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
      mat_acc45_post_rfx <- megamatrix_delay(F_params=F_params,
                                             M_params=M_params,
                                             twosex=F,
                                             grow_perturb=0,
                                             surv_perturb=0,
                                             flow_perturb=0,
                                             fert_perturb=0,
                                             viab_perturb=0,
                                             zpptgrow=climate_acc_45_values[l,5],
                                             ztempgrow=climate_acc_45_values[l,6],
                                             zpptdorm=climate_acc_45_values[l,7],
                                             ztempdorm=climate_acc_45_values[l,8],
                                             rfx=rfx_fun())$MEGAmat
      popbio::lambda(mat_acc45_post_rfx)
      
    }
)

geo_lambba_acc45_fd<-data.frame(lambda_acc45)

lambda_post_acc85<-matrix(NA,nrow=n_post_draws,ncol=nrow(climate_acc_85_values))
dim(lambda_post_acc85)
# lambda_dorm_post<-c()

#### ACCESS 85----
F_params <- M_params <- list()
Timeacc85<-system.time(
  lambda_acc85<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(climate_acc_85_values),.combine='c') %dopar% {
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
      mat_acc85_post_rfx <- megamatrix_delay(F_params=F_params,
                                             M_params=M_params,
                                             twosex=F,
                                             grow_perturb=0,
                                             surv_perturb=0,
                                             flow_perturb=0,
                                             fert_perturb=0,
                                             viab_perturb=0,
                                             zpptgrow=climate_acc_85_values[l,5],
                                             ztempgrow=climate_acc_85_values[l,6],
                                             zpptdorm=climate_acc_85_values[l,7],
                                             ztempdorm=climate_acc_85_values[l,8],
                                             rfx=rfx_fun())$MEGAmat
      popbio::lambda(mat_acc85_post_rfx)
      
    }
)

geolambba_acc85_fd<-data.frame(lambda_acc85)
# write.csv(data_acc85_fd,"C:/Users/jm200/Documents/Projection/2024/Geoprojectfd/data_acc85_fd.csv")
# write_rds(data_acc85_fd,"C:/Users/jm200/Documents/Projection/2024/Geoprojectfd/data_acc85_fd.rds")

#### CESM 45----
lambda_post_ces45<-matrix(NA,nrow=n_post_draws,ncol=nrow(climate_ces_45_values))
F_params <- M_params <- list()
Timeces45<-system.time(
  lambda_ces45<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(climate_ces_45_values),.combine='c') %dopar% {
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
      mat_ces45_post_rfx <- megamatrix_delay(F_params=F_params,
                                             M_params=M_params,
                                             twosex=F,
                                             grow_perturb=0,
                                             surv_perturb=0,
                                             flow_perturb=0,
                                             fert_perturb=0,
                                             viab_perturb=0,
                                             zpptgrow=climate_ces_45_values[l,5],
                                             ztempgrow=climate_ces_45_values[l,6],
                                             zpptdorm=climate_ces_45_values[l,7],
                                             ztempdorm=climate_ces_45_values[l,8],
                                             rfx=rfx_fun())$MEGAmat
      popbio::lambda(mat_ces45_post_rfx)
      
    }
)

geo_lambba_ces45_fd<-data.frame(lambda_ces45)

lambda_post_ces85<-matrix(NA,nrow=n_post_draws,ncol=nrow(climate_ces_85_values))
dim(lambda_post_ces85)
# lambda_dorm_post<-c()

#### CESM 85-----
F_params <- M_params <- list()
Timeces85<-system.time(
  lambda_ces85<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(climate_ces_85_values),.combine='c') %dopar% {
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
      mat_ces85_post_rfx <- megamatrix_delay(F_params=F_params,
                                             M_params=M_params,
                                             twosex=F,
                                             grow_perturb=0,
                                             surv_perturb=0,
                                             flow_perturb=0,
                                             fert_perturb=0,
                                             viab_perturb=0,
                                             zpptgrow=climate_ces_85_values[l,5],
                                             ztempgrow=climate_ces_85_values[l,6],
                                             zpptdorm=climate_ces_85_values[l,7],
                                             ztempdorm=climate_ces_85_values[l,8],
                                             rfx=rfx_fun())$MEGAmat
      popbio::lambda(mat_ces85_post_rfx)
      
    }
)

geolambba_ces85_fd<-data.frame(lambda_ces85)

#### CMCM 45----
lambda_post_cmc45<-matrix(NA,nrow=n_post_draws,ncol=nrow(climate_cmc_45_values))
#dim(lambda_post_cmc5)
F_params <- M_params <- list()
Timecmc45<-system.time(
  lambda_cmc45<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(climate_cmc_45_values),.combine='c') %dopar% {
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
      mat_cmc45_post_rfx <- megamatrix_delay(F_params=F_params,
                                             M_params=M_params,
                                             twosex=F,
                                             grow_perturb=0,
                                             surv_perturb=0,
                                             flow_perturb=0,
                                             fert_perturb=0,
                                             viab_perturb=0,
                                             zpptgrow=climate_cmc_45_values[l,5],
                                             ztempgrow=climate_cmc_45_values[l,6],
                                             zpptdorm=climate_cmc_45_values[l,7],
                                             ztempdorm=climate_cmc_45_values[l,8],
                                             rfx=rfx_fun())$MEGAmat
      popbio::lambda(mat_cmc45_post_rfx)
      
    }
)

geo_lambba_cmc45_fd<-data.frame(lambda_cmc45)

#### CMCM 85----
lambda_post_cmc85<-matrix(NA,nrow=n_post_draws,ncol=nrow(climate_cmc_85_values))
#dim(lambda_post_cmc5)
F_params <- M_params <- list()
Timecmc5<-system.time(
  lambda_cmc85<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(climate_cmc_85_values),.combine='c') %dopar% {
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
      mat_cmc85_post_rfx <- megamatrix_delay(F_params=F_params,
                                             M_params=M_params,
                                             twosex=F,
                                             grow_perturb=0,
                                             surv_perturb=0,
                                             flow_perturb=0,
                                             fert_perturb=0,
                                             viab_perturb=0,
                                             zpptgrow=climate_cmc_85_values[l,5],
                                             ztempgrow=climate_cmc_85_values[l,6],
                                             zpptdorm=climate_cmc_85_values[l,7],
                                             ztempdorm=climate_cmc_85_values[l,8],
                                             rfx=rfx_fun())$MEGAmat
      popbio::lambda(mat_cmc85_post_rfx)
      
    }
)

geolambba_cmc85_fd<-data.frame(lambda_cmc85)

#### MIROC 45----
lambda_post_miroc45<-matrix(NA,nrow=n_post_draws,ncol=nrow(climate_cmc_45_values))
#dim(lambda_post_cmc5)
F_params <- M_params <- list()
Timemiroc45<-system.time(
  lambda_miroc45<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(climate_miroc45_values),.combine='c') %dopar% {
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
      mat_miroc45_post_rfx <- megamatrix_delay(F_params=F_params,
                                             M_params=M_params,
                                             twosex=F,
                                             grow_perturb=0,
                                             surv_perturb=0,
                                             flow_perturb=0,
                                             fert_perturb=0,
                                             viab_perturb=0,
                                             zpptgrow=climate_miroc45_values[l,5],
                                             ztempgrow=climate_miroc45_values[l,6],
                                             zpptdorm=climate_miroc45_values[l,7],
                                             ztempdorm=climate_miroc45_values[l,8],
                                             rfx=rfx_fun())$MEGAmat
      popbio::lambda(mat_miroc45_post_rfx)
      
    }
)

geo_lambda_miroc45_fd<-data.frame(lambda_cmc45)

#### MIROC 85----
lambda_post_miroc85<-matrix(NA,nrow=n_post_draws,ncol=nrow(climate_miroc_85_values))
#dim(lambda_post_cmc5)
F_params <- M_params <- list()
Timemiroc85<-system.time(
  lambda_miroc85<-foreach(p=1:n_post_draws,.combine='rbind') %:% 
    foreach(l=1:nrow(climate_miroc85_values),.combine='c') %dopar% {
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
      mat_miroc85_post_rfx <- megamatrix_delay(F_params=F_params,
                                             M_params=M_params,
                                             twosex=F,
                                             grow_perturb=0,
                                             surv_perturb=0,
                                             flow_perturb=0,
                                             fert_perturb=0,
                                             viab_perturb=0,
                                             zpptgrow=climate_miroc85_values[l,5],
                                             ztempgrow=climate_miroc85_values[l,6],
                                             zpptdorm=climate_miroc85_values[l,7],
                                             ztempdorm=climate_miroc85_values[l,8],
                                             rfx=rfx_fun())$MEGAmat
      popbio::lambda(mat_miroc85_post_rfx)
      
    }
)

geo_lambba_miroc85_fd<-data.frame(lambda_miroc85)

# Fig 4. Geography ----
## Download occurrence data from gbif for Poa arachnifera -----
# dir.create("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/data/occurence")
poar_occ_raw <- dismo::gbif(genus="Poa",species="arachnifera",download=TRUE) 
head(poar_occ_raw) 
poar_occ <- subset(poar_occ_raw,(!is.na(lat))&(!is.na(lon))) # here we remove erroneous coordinates, where either the latitude or longitude is missing
cat(nrow(poar_occ_raw)-nrow(poar_occ), "records are removed") # Show the number of records that are removed from the dataset. 
poar_occ %>% 
  dplyr::select(country,lon, lat,year)%>% 
  dplyr::rename(Longitude=lon,Latitude=lat) %>% 
  filter(year %in% (1901:2024) & as.numeric(Longitude >=-102) &  as.numeric(Longitude <=-95) & country=="United States") %>% 
  unique() %>% 
  arrange(Latitude)->gbif
# coordinates(gbif) <- ~ Longitude + Latitude
# CRS1 <- CRS("+init=epsg:4326") # WGS 84
# crs(gbif) <- CRS1

## Read the lambda for all iterations
### Two sex models
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
geo_lambbacurrent_fd <- 
geo_lambbapast_fd <- 
geo_lambda_miroc45_fd<-
geo_lambba_miroc85_fd <- 
geo_lambba_cmc45_fd <- 
geolambba_cmc85_fd <- 
geo_lambba_acc45_fd <- 
geolambba_acc85_fd <- 
geo_lambba_ces45_fd <- 
geolambba_ces85_fd <- 

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

### Present conditions -----
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

### Future conditions -----
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
# Peprare the data for the plot
geo_prlambda_past_fd_twosex<-left_join(lam_prob_past_fd,geo_prlambda_past,by=c("x","y"))
geo_prlambda_past_fd_twosex$Prlambda<-geo_prlambda_past_fd_twosex$Prlambda.x-geo_prlambda_past_fd_twosex$Prlambda.y
geo_prlambda_current_fd_twosex<-left_join(lam_prob_current_fd,geo_prlambda_current,by=c("x","y"))
geo_prlambda_current_fd_twosex$Prlambda<-geo_prlambda_current_fd_twosex$Prlambda.x-geo_prlambda_current_fd_twosex$Prlambda.y
geo_prlambda_miroc45_fd_twosex<-left_join(lam_prob_miroc45_fd,geo_prlambda_miroc45,by=c("x","y"))
geo_prlambda_miroc45_fd_twosex$Prlambda<-geo_prlambda_miroc45_fd_twosex$Prlambda.x-geo_prlambda_miroc45_fd_twosex$Prlambda.y
geo_prlambda_miroc85_fd_twosex<-left_join(lam_prob_miroc85_fd,geo_prlambda_miroc85,by=c("x","y"))
geo_prlambda_miroc85_fd_twosex$Prlambda<-geo_prlambda_miroc85_fd_twosex$Prlambda.x-geo_prlambda_miroc85_fd_twosex$Prlambda.y

dat_miroc_map<-rbind(geo_prlambda_past[,c("x","y","Prlambda")],lam_prob_past_fd[,c("x","y","Prlambda")],geo_prlambda_past_fd_twosex[,c("x","y","Prlambda")],geo_prlambda_current[,c("x","y","Prlambda")],lam_prob_current_fd[,c("x","y","Prlambda")],geo_prlambda_current_fd_twosex[,c("x","y","Prlambda")],geo_prlambda_miroc45[,c("x","y","Prlambda")],lam_prob_miroc45_fd[,c("x","y","Prlambda")],geo_prlambda_miroc45_fd_twosex[,c("x","y","Prlambda")],geo_prlambda_miroc85[,c("x","y","Prlambda")],lam_prob_miroc85_fd[,c("x","y","Prlambda")],geo_prlambda_miroc85_fd_twosex[,c("x","y","Prlambda")])

# dim(dat_miroc_map)
dat_miroc_map$Time<-c(rep("Past",nrow(geo_prlambda_past)+nrow(lam_prob_past_fd)+nrow(geo_prlambda_past_fd_twosex)),rep("Present",3*nrow(geo_prlambda_current)),rep("RCP45",3*nrow(geo_prlambda_miroc45)),rep("RCP85",3*nrow(geo_prlambda_miroc85)))
dat_miroc_map$model<-c(rep("B",nrow(geo_prlambda_past)),rep("A",nrow(lam_prob_past_fd)),rep("C",nrow(geo_prlambda_past_fd_twosex)),
                       rep("B",nrow(geo_prlambda_current)),rep("A",nrow(lam_prob_current_fd)),rep("C",nrow(geo_prlambda_current_fd_twosex)),
                       rep("B",nrow(geo_prlambda_miroc45)),rep("A",nrow(lam_prob_miroc45_fd)),rep("C",nrow(geo_prlambda_miroc45_fd_twosex)),
                       
                       rep("B",nrow(geo_prlambda_miroc85)),rep("A",nrow(lam_prob_miroc85_fd)),rep("C",nrow(geo_prlambda_miroc85_fd_twosex)) 
)

model_names <- c("A"="Female dominant model (F)",
                 "B"="Two-sex model (FM)",
                 "C"="F - FM"
)
Time_names <- c("Past"="Past",
                "Present"="Present",
                "RCP45"="RCP 4.5",
                "RCP85"="RCP 8.5"
)

occplot<-data.frame(gbif,Time=c(rep("Past",3*nrow(gbif)),rep("Present",3*nrow(gbif)),rep("RCP45",3*nrow(gbif)),rep("RCP85",3*nrow(gbif))))
occplot$model<-c(rep("A",nrow(gbif)),rep("B",nrow(gbif)),rep("C",nrow(gbif)))

                    
Fig_geoPrlambdamiroc<-ggplot()+
  geom_tile(data = subset(dat_miroc_map, model == "A"), aes(x = x, y = y, fill = Prlambda))+
  #geom_point(data = subset(occplot, Time_names=="Present" & model == "A"), aes(x = Longitude, y = Latitude),size=0.01)+
  # facet_grid(model~Time,scales='free_x', space='free_x',labeller = labeller(model = model_names,Time=Time_names)) +
  labs(x = "Longitude", y = "Latitude")+
  scale_fill_gradientn(
    name = expression("Pr " (lambda[F]) > 1),
    colours = terrain.colors(100),
    na.value = "transparent",
    breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
    limits=c(0,1))+
  new_scale_fill() +
  geom_tile(data = subset(dat_miroc_map, model == "B"), aes(x = x, y = y, fill = Prlambda))+
  #geom_point(data = subset(occplot, Time_names=="Present" & model == "B"), aes(x = Longitude, y = Latitude),size=0.01)+
  # facet_grid(model~Time,scales='free_x', space='free_x',labeller = labeller(model = model_names,Time=Time_names)) +
  labs(x = "Longitude", y = "Latitude")+
  scale_fill_gradientn(
    name = expression("Pr " (lambda[FM]) > 1),
    colours = terrain.colors(100),
    na.value = "transparent",
    breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
    limits=c(0,1))+
  new_scale_fill() +
  geom_tile(data = subset(dat_miroc_map, model == "C"), aes(x = x, y = y, fill = Prlambda))+
  # geom_point(data = occplot, aes(x = Longitude, y = Latitude),Time="Present",size=0.1)+
  scale_fill_gradientn(
    name = expression(paste(Delta,"Pr " (lambda[F-FM]))> 1),
    colours = colorspace::diverge_hcl(7),
    na.value = "transparent",
    breaks = seq(-0.5, 0.5, length.out=5),labels=c(-0.50,-0.25,0.00,0.25,0.50),
    limits=c(-0.5,0.5))+
  facet_grid(model~Time,labeller = labeller(model = model_names,Time=Time_names))+
  geom_point(data = subset(occplot, Time=="Present" & model==c("A","B")), aes(x = Longitude, y = Latitude),size=1,alpha = 1/5)+
  theme_light()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 7),
        strip.text.x = element_text(
          size = 10, color = "grey40", face = "bold"
        ),
        strip.text.y = element_text(
          size = 10, color = "grey40", face = "bold"),
        strip.background = element_rect(
          color="black", fill="white", size=0.5, linetype="solid"
        )
  )


ggsave("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Fig_geoPrlambda_miroc.pdf", Fig_geoPrlambdamiroc, width =9, height = 10)

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

geo_prlambda_acc45_fd_twosex<-left_join(lam_prob_acc45_fd,geo_prlambdaacc_45,by=c("x","y"))
geo_prlambda_acc45_fd_twosex$Prlambda<-geo_prlambda_acc45_fd_twosex$Prlambda.x-geo_prlambda_acc45_fd_twosex$Prlambda.y
geo_prlambda_acc85_fd_twosex<-left_join(lam_prob_acc85_fd,geo_prlambdaacc_85,by=c("x","y"))
geo_prlambda_acc85_fd_twosex$Prlambda<-geo_prlambda_acc85_fd_twosex$Prlambda.x-geo_prlambda_acc85_fd_twosex$Prlambda.y

dat_acc_map<-rbind(geo_prlambda_past[,c("x","y","Prlambda")],lam_prob_past_fd[,c("x","y","Prlambda")],geo_prlambda_past_fd_twosex[,c("x","y","Prlambda")],geo_prlambda_current[,c("x","y","Prlambda")],lam_prob_current_fd[,c("x","y","Prlambda")],geo_prlambda_current_fd_twosex[,c("x","y","Prlambda")],geo_prlambdaacc_45[,c("x","y","Prlambda")],lam_prob_acc45_fd[,c("x","y","Prlambda")],geo_prlambda_acc45_fd_twosex[,c("x","y","Prlambda")],geo_prlambdaacc_85[,c("x","y","Prlambda")],lam_prob_acc85_fd[,c("x","y","Prlambda")],geo_prlambda_acc85_fd_twosex[,c("x","y","Prlambda")])

dat_acc_map$Time<-c(rep("Past",nrow(geo_prlambda_past)+nrow(lam_prob_past_fd)+nrow(geo_prlambda_past_fd_twosex)),rep("Present",3*nrow(geo_prlambda_current)),rep("RCP45",3*nrow(geo_prlambdaacc_45)),rep("RCP85",3*nrow(geo_prlambdaacc_85)))

dat_acc_map$model<-c(rep("B",nrow(geo_prlambda_past)),rep("A",nrow(lam_prob_past_fd)),rep("C",nrow(geo_prlambda_past_fd_twosex)),
                     rep("B",nrow(geo_prlambda_current)),rep("A",nrow(lam_prob_current_fd)),rep("C",nrow(geo_prlambda_current_fd_twosex)),
                     rep("B",nrow(geo_prlambdaacc_45)),rep("A",nrow(lam_prob_acc45_fd)),rep("C",nrow(geo_prlambda_acc45_fd_twosex)),
                     
                     rep("B",nrow(geo_prlambdaacc_85)),rep("A",nrow(lam_prob_acc85_fd)),rep("C",nrow(geo_prlambda_acc85_fd_twosex)) 
)

Fig_geoPrlambdaacc<-ggplot()+
  geom_tile(data = subset(dat_acc_map, model == "A"), aes(x = x, y = y, fill = Prlambda))+
  # facet_grid(model~Time,scales='free_x', space='free_x',labeller = labeller(model = model_names,Time=Time_names)) +
  labs(x = "Longitude", y = "Latitude")+
  scale_fill_gradientn(
    name = expression("Pr " (lambda[F]) > 1),
    colours = terrain.colors(100),
    na.value = "transparent",
    breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
    limits=c(0,1))+
  new_scale_fill() +
  geom_tile(data = subset(dat_acc_map, model == "B"), aes(x = x, y = y, fill = Prlambda))+
  # facet_grid(model~Time,scales='free_x', space='free_x',labeller = labeller(model = model_names,Time=Time_names)) +
  labs(x = "Longitude", y = "Latitude")+
  scale_fill_gradientn(
    name = expression("Pr " (lambda[FM]) > 1),
    colours = terrain.colors(100),
    na.value = "transparent",
    breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
    limits=c(0,1))+
  new_scale_fill() +
  geom_tile(data = subset(dat_acc_map, model == "C"), aes(x = x, y = y, fill = Prlambda))+
  scale_fill_gradientn(
    name = expression(paste(Delta,"Pr " (lambda[F-FM]))> 1),
    colours = colorspace::diverge_hcl(7),
    na.value = "transparent",
    breaks = seq(-0.5, 0.5, length.out=5),labels=c(-0.5,-0.25,0.00,0.25,0.5),
    limits=c(-0.5,0.5))+
  facet_grid(model~Time,labeller = labeller(model = model_names,Time=Time_names))+
  geom_point(data = subset(occplot, Time=="Present" & model==c("A","B")), aes(x = Longitude, y = Latitude),size=1.5,alpha = 1/5)+
  theme_light()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 7),
        strip.text.x = element_text(
          size = 10, color = "grey40", face = "bold"
        ),
        strip.text.y = element_text(
          size = 10, color = "grey40", face = "bold"),
        strip.background = element_rect(
          color="black", fill="white", size=0.5, linetype="solid"
        )
  )

ggsave("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Fig_geoPrlambdaacc.pdf", Fig_geoPrlambdaacc, width = 9, height = 10)

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

geo_prlambda_cmc45_fd_twosex<-left_join(lam_prob_cmc45_fd,geo_prlambdacmc_45,by=c("x","y"))
geo_prlambda_cmc45_fd_twosex$Prlambda<-geo_prlambda_cmc45_fd_twosex$Prlambda.x-geo_prlambda_cmc45_fd_twosex$Prlambda.y

geo_prlambda_cmc85_fd_twosex<-left_join(lam_prob_cmc85_fd,geo_prlambdacmc_85,by=c("x","y"))
geo_prlambda_cmc85_fd_twosex$Prlambda<-geo_prlambda_cmc85_fd_twosex$Prlambda.x-geo_prlambda_cmc85_fd_twosex$Prlambda.y

dat_cmc_map<-rbind(geo_prlambda_past[,c("x","y","Prlambda")],lam_prob_past_fd[,c("x","y","Prlambda")],geo_prlambda_past_fd_twosex[,c("x","y","Prlambda")],geo_prlambda_current[,c("x","y","Prlambda")],lam_prob_current_fd[,c("x","y","Prlambda")],geo_prlambda_current_fd_twosex[,c("x","y","Prlambda")],geo_prlambdacmc_45[,c("x","y","Prlambda")],lam_prob_cmc45_fd[,c("x","y","Prlambda")],geo_prlambda_cmc45_fd_twosex[,c("x","y","Prlambda")],geo_prlambdacmc_85[,c("x","y","Prlambda")],lam_prob_cmc85_fd[,c("x","y","Prlambda")],geo_prlambda_cmc85_fd_twosex[,c("x","y","Prlambda")])



dat_cmc_map$Time<-c(rep("Past",nrow(geo_prlambda_past)+nrow(lam_prob_past_fd)+nrow(geo_prlambda_past_fd_twosex)),rep("Present",3*nrow(geo_prlambda_current)),rep("RCP45",3*nrow(geo_prlambdacmc_45)),rep("RCP85",3*nrow(geo_prlambdacmc_85)))

dat_cmc_map$model<-c(rep("B",nrow(geo_prlambda_past)),rep("A",nrow(lam_prob_past_fd)),rep("C",nrow(geo_prlambda_past_fd_twosex)),
                     rep("B",nrow(geo_prlambda_current)),rep("A",nrow(lam_prob_current_fd)),rep("C",nrow(geo_prlambda_current_fd_twosex)),
                     rep("B",nrow(geo_prlambdacmc_45)),rep("A",nrow(lam_prob_cmc45_fd)),rep("C",nrow(geo_prlambda_cmc45_fd_twosex)),
                     
                     rep("B",nrow(geo_prlambdacmc_85)),rep("A",nrow(lam_prob_cmc85_fd)),rep("C",nrow(geo_prlambda_cmc85_fd_twosex)) 
)


Fig_geoPrlambdacmc<-ggplot()+
  geom_tile(data = subset(dat_cmc_map, model == "A"), aes(x = x, y = y, fill = Prlambda))+
  # facet_grid(model~Time,scales='free_x', space='free_x',labeller = labeller(model = model_names,Time=Time_names)) +
  labs(x = "Longitude", y = "Latitude")+
  scale_fill_gradientn(
    name = expression("Pr " (lambda[F]) > 1),
    colours = terrain.colors(100),
    na.value = "transparent",
    breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
    limits=c(0,1))+
  new_scale_fill() +
  geom_tile(data = subset(dat_cmc_map, model == "B"), aes(x = x, y = y, fill = Prlambda))+
  # facet_grid(model~Time,scales='free_x', space='free_x',labeller = labeller(model = model_names,Time=Time_names)) +
  labs(x = "Longitude", y = "Latitude")+
  scale_fill_gradientn(
    name = expression("Pr " (lambda[FM]) > 1),
    colours = terrain.colors(100),
    na.value = "transparent",
    breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
    limits=c(0,1))+
  new_scale_fill() +
  geom_tile(data = subset(dat_cmc_map, model == "C"), aes(x = x, y = y, fill = Prlambda))+
  scale_fill_gradientn(
    name = expression(paste(Delta,"Pr " (lambda[F-FM]))> 1),
    colours = colorspace::diverge_hcl(7),
    na.value = "transparent",
    breaks = seq(-0.5, 0.5, length.out=5),labels=c(-0.5,-0.25,0.00,0.25,0.5),
    limits=c(-0.5,0.5))+
  facet_grid(model~Time,labeller = labeller(model = model_names,Time=Time_names))+
  geom_point(data = subset(occplot, Time=="Present" & model==c("A","B")), aes(x = Longitude, y = Latitude),size=1.5,alpha = 1/5)+
  theme_light()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 7),
        strip.text.x = element_text(
          size = 12, color = "grey40", face = "bold"
        ),
        strip.text.y = element_text(
          size = 12, color = "grey40", face = "bold"),
        strip.background = element_rect(
          color="black", fill="white", size=0.5, linetype="solid"
        )
  )

ggsave("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Fig_geoPrlambdacmc.pdf", Fig_geoPrlambdacmc, width = 9, height = 10)

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

geo_prlambda_ces45_fd_twosex<-left_join(lam_prob_ces45_fd,geo_prlambdaces_45,by=c("x","y"))
geo_prlambda_ces45_fd_twosex$Prlambda<-geo_prlambda_ces45_fd_twosex$Prlambda.x-geo_prlambda_ces45_fd_twosex$Prlambda.y

geo_prlambda_ces85_fd_twosex<-left_join(lam_prob_ces85_fd,geo_prlambdaces_85,by=c("x","y"))
geo_prlambda_ces85_fd_twosex$Prlambda<-geo_prlambda_ces85_fd_twosex$Prlambda.x-geo_prlambda_ces85_fd_twosex$Prlambda.y

dat_ces_map<-rbind(geo_prlambda_past[,c("x","y","Prlambda")],lam_prob_past_fd[,c("x","y","Prlambda")],geo_prlambda_past_fd_twosex[,c("x","y","Prlambda")],geo_prlambda_current[,c("x","y","Prlambda")],lam_prob_current_fd[,c("x","y","Prlambda")],geo_prlambda_current_fd_twosex[,c("x","y","Prlambda")],geo_prlambdaces_45[,c("x","y","Prlambda")],lam_prob_ces45_fd[,c("x","y","Prlambda")],geo_prlambda_ces45_fd_twosex[,c("x","y","Prlambda")],geo_prlambdaces_85[,c("x","y","Prlambda")],lam_prob_ces85_fd[,c("x","y","Prlambda")],geo_prlambda_ces85_fd_twosex[,c("x","y","Prlambda")])

dim(dat_ces_map)

dat_ces_map$Time<-c(rep("Past",nrow(geo_prlambda_past)+nrow(lam_prob_past_fd)+nrow(geo_prlambda_past_fd_twosex)),rep("Present",3*nrow(geo_prlambda_current)),rep("RCP45",3*nrow(geo_prlambdaces_45)),rep("RCP85",3*nrow(geo_prlambdaces_85)))

dat_ces_map$model<-c(rep("B",nrow(geo_prlambda_past)),rep("A",nrow(lam_prob_past_fd)),rep("C",nrow(geo_prlambda_past_fd_twosex)),
                     rep("B",nrow(geo_prlambda_current)),rep("A",nrow(lam_prob_current_fd)),rep("C",nrow(geo_prlambda_current_fd_twosex)),
                     rep("B",nrow(geo_prlambdaces_45)),rep("A",nrow(lam_prob_ces45_fd)),rep("C",nrow(geo_prlambda_ces45_fd_twosex)),
                     
                     rep("B",nrow(geo_prlambdaces_85)),rep("A",nrow(lam_prob_ces85_fd)),rep("C",nrow(geo_prlambda_ces85_fd_twosex)) 
)

Fig_geoPrlambdaces<-ggplot()+
  geom_tile(data = subset(dat_ces_map, model == "A"), aes(x = x, y = y, fill = Prlambda))+
  # facet_grid(model~Time,scales='free_x', space='free_x',labeller = labeller(model = model_names,Time=Time_names)) +
  labs(x = "Longitude", y = "Latitude")+
  scale_fill_gradientn(
    name = expression("Pr " (lambda[F]) > 1),
    colours = terrain.colors(100),
    na.value = "transparent",
    breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
    limits=c(0,1))+
  new_scale_fill() +
  geom_tile(data = subset(dat_ces_map, model == "B"), aes(x = x, y = y, fill = Prlambda))+
  # facet_grid(model~Time,scales='free_x', space='free_x',labeller = labeller(model = model_names,Time=Time_names)) +
  labs(x = "Longitude", y = "Latitude")+
  scale_fill_gradientn(
    name = expression("Pr " (lambda[FM]) > 1),
    colours = terrain.colors(100),
    na.value = "transparent",
    breaks = seq(0, 1, length.out=5),labels=c(0.00, 0.25, 0.50, 0.75, 1.00),
    limits=c(0,1))+
  new_scale_fill() +
  geom_tile(data = subset(dat_ces_map, model == "C"), aes(x = x, y = y, fill = Prlambda))+
  scale_fill_gradientn(
    name = expression(paste(Delta,"Pr " (lambda[F-FM]))> 1),
    colours = colorspace::diverge_hcl(7),
    na.value = "transparent",
    breaks = seq(-0.5, 0.5, length.out=5),labels=c(-0.5,-0.25,0.00,0.25,0.5),
    limits=c(-0.5,0.5))+
  facet_grid(model~Time,labeller = labeller(model = model_names,Time=Time_names))+
  geom_point(data = subset(occplot, Time=="Present" & model==c("A","B")), aes(x = Longitude, y = Latitude),size=1.5,alpha = 1/5)+
  theme_light()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 7),
        strip.text.x = element_text(
          size = 10, color = "grey40", face = "bold"
        ),
        strip.text.y = element_text(
          size = 10, color = "grey40", face = "bold"),
        strip.background = element_rect(
          color="black", fill="white", size=0.5, linetype="solid"
        )
  )

ggsave("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Fig_geoPrlambda_ces.pdf", Fig_geoPrlambdaces,  width = 9, height = 10)

# Figure: density plot -----
d_past_lambda_fd <- density(unique(lam_prob_past_fd$Prlambda),na.rm=TRUE)
d_past_lambda_2sex <- density(unique(geo_prlambda_past$Prlambda),na.rm=TRUE)
d_past_lambda_fd_2sex <- density(unique(geo_prlambda_past_fd_twosex$Prlambda),na.rm=TRUE)

d_current_lambda_fd <- density(unique(lam_prob_current_fd$Prlambda),na.rm=TRUE)
d_current_lambda_2sex <- density(unique(geo_prlambda_current$Prlambda),na.rm=TRUE)
d_current_lambda_fd_2sex <- density(unique(geo_prlambda_current_fd_twosex$Prlambda),na.rm=TRUE)

d_miroc45_lambda_fd <- density(unique(lam_prob_miroc45_fd$Prlambda),na.rm=TRUE)
d_miroc45_lambda_2sex <- density(unique(geo_prlambda_miroc45$Prlambda),na.rm=TRUE)
d_miroc45_lambda_fd_2sex <- density(unique(geo_prlambda_miroc45_fd_twosex$Prlambda),na.rm=TRUE)
d_miroc85_lambda_fd <- density(unique(lam_prob_miroc85_fd$Prlambda),na.rm=TRUE)
d_miroc85_lambda_2sex <- density(unique(geo_prlambda_miroc85$Prlambda),na.rm=TRUE)
d_miroc85_lambda_fd_2sex <- density(unique(geo_prlambda_miroc85_fd_twosex$Prlambda),na.rm=TRUE)

d_acc45_lambda_fd <- density(unique(lam_prob_acc45_fd$Prlambda),na.rm=TRUE)
d_acc45_lambda_2sex <- density(unique(geo_prlambdaacc_45$Prlambda),na.rm=TRUE)
d_acc45_lambda_fd_2sex <- density(unique(geo_prlambda_acc45_fd_twosex$Prlambda),na.rm=TRUE)
d_acc85_lambda_fd <- density(unique(lam_prob_acc85_fd$Prlambda),na.rm=TRUE)
d_acc85_lambda_2sex <- density(unique(geo_prlambdaacc_85$Prlambda),na.rm=TRUE)
d_acc85_lambda_fd_2sex <- density(unique(geo_prlambda_acc85_fd_twosex$Prlambda),na.rm=TRUE)


d_cmc45_lambda_fd <- density(unique(lam_prob_cmc45_fd$Prlambda),na.rm=TRUE)
d_cmc45_lambda_2sex <- density(unique(geo_prlambdacmc_45$Prlambda),na.rm=TRUE)
d_cmc45_lambda_fd_2sex <- density(unique(geo_prlambda_cmc45_fd_twosex$Prlambda),na.rm=TRUE)
d_cmc85_lambda_fd <- density(unique(lam_prob_cmc85_fd$Prlambda),na.rm=TRUE)
d_cmc85_lambda_2sex <- density(unique(geo_prlambdacmc_85$Prlambda),na.rm=TRUE)
d_cmc85_lambda_fd_2sex <- density(unique(geo_prlambda_cmc85_fd_twosex$Prlambda),na.rm=TRUE)

d_ces45_lambda_fd <- density(unique(lam_prob_ces45_fd$Prlambda),na.rm=TRUE)
d_ces45_lambda_2sex <- density(unique(geo_prlambdaces_45$Prlambda),na.rm=TRUE)
d_ces45_lambda_fd_2sex <- density(unique(geo_prlambda_ces45_fd_twosex$Prlambda),na.rm=TRUE)
d_ces85_lambda_fd <- density(unique(lam_prob_ces85_fd$Prlambda),na.rm=TRUE)
d_ces85_lambda_2sex <- density(unique(geo_prlambdaces_85$Prlambda),na.rm=TRUE)
d_ces85_lambda_fd_2sex <- density(unique(geo_prlambda_ces85_fd_twosex$Prlambda),na.rm=TRUE)



pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Densityplot_lambda_Miroc5.pdf",width=4,height=8,useDingbats = F)
par(mar=c(5,5,1.5,0.5),mfrow=c(2,1))
plot(d_past_lambda_fd, lwd = 2, main = "", xlab =expression("Pr " (lambda[F]) > 1),ylim=c(0,3),xlim=c(0,1),col = "#7570B3")
# mtext( "A",side = 3, adj = 0,cex=1.25)
abline(v=mean(lam_prob_past_fd$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_lambda_fd, lwd = 2,col = "#0072B2")
abline(v=mean(lam_prob_current_fd$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_miroc45_lambda_fd, lwd = 2,col = "#E6AB02")
abline(v=mean(lam_prob_miroc45_fd$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_miroc85_lambda_fd, lwd = 2,col = "#E7298A")
abline(v=mean(lam_prob_miroc85_fd$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")


legend(0.01,3,bty="n",legend=c("Past","Present","RCP 4.5","RCP 8.5"),lwd=c(2,2,2,2),lty=1,col=c("#7570B3","#0072B2","#E6AB02","#E7298A"),cex=0.8)

plot(d_past_lambda_2sex, lwd = 2, main = "", xlab = expression("Pr " (lambda[FM]) > 1),,xlim=c(0,1),ylim=c(0,3),col = "#7570B3")
# mtext( "A",side = 3, adj = 0,cex=1.25)
abline(v=mean(geo_prlambda_past$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_lambda_2sex, lwd = 2,col = "#0072B2")
abline(v=mean(geo_prlambda_current$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_miroc45_lambda_2sex, lwd = 2,col = "#E6AB02")
abline(v=mean(geo_prlambda_miroc45$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_miroc85_lambda_2sex, lwd = 2,col = "#E7298A")
abline(v=mean(geo_prlambda_miroc85$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")

dev.off()

pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Densityplot_lambda_acc.pdf",width=4,height=8,useDingbats = F)
par(mar=c(5,5,1.5,0.5),mfrow=c(2,1))
plot(d_past_lambda_fd, lwd = 2, main = "", xlab = expression("Pr " (lambda[F]) > 1),xlim=c(0,1),ylim=c(0,2),col = "#7570B3")
# mtext( "A",side = 3, adj = 0,cex=1.25)
abline(v=mean(lam_prob_past_fd$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_lambda_fd, lwd = 2,col = "#0072B2")
abline(v=mean(lam_prob_current_fd$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_acc45_lambda_fd, lwd = 2,col = "#E6AB02")
abline(v=mean(lam_prob_acc45_fd$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_acc85_lambda_fd, lwd = 2,col = "#E7298A")
abline(v=mean(lam_prob_acc85_fd$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")


legend(0.01,2.2,bty="n",legend=c("Past","Present","RCP 4.5","RCP 8.5"),lwd=c(2,2,2,2),lty=1,col=c("#7570B3","#0072B2","#E6AB02","#E7298A"),cex=0.8)

plot(d_past_lambda_2sex, lwd = 2, main = "", xlab = expression("Pr " (lambda[FM]) > 1),,xlim=c(0,1),ylim=c(0,3),col = "#7570B3")
# mtext( "A",side = 3, adj = 0,cex=1.25)
abline(v=mean(geo_prlambda_past$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_lambda_2sex, lwd = 2,col = "#0072B2")
abline(v=mean(geo_prlambda_current$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_acc45_lambda_2sex, lwd = 2,col = "#E6AB02")
abline(v=mean(geo_prlambdaacc_45$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_acc85_lambda_2sex, lwd = 2,col = "#E7298A")
abline(v=mean(geo_prlambdaacc_85$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")

dev.off()




pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Densityplot_lambda_cmc.pdf",width=4,height=8,useDingbats = F)
par(mar=c(5,5,1.5,0.5),mfrow=c(2,1))
plot(d_past_lambda_fd, lwd = 2, main = "", xlab = expression("Pr " (lambda[F]) > 1),xlim=c(0,1),ylim=c(0,2),col = "#7570B3")
# mtext( "A",side = 3, adj = 0,cex=1.25)
abline(v=mean(lam_prob_past_fd$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_lambda_fd, lwd = 2,col = "#0072B2")
abline(v=mean(lam_prob_current_fd$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_cmc45_lambda_fd, lwd = 2,col = "#E6AB02")
abline(v=mean(lam_prob_cmc45_fd$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_cmc85_lambda_fd, lwd = 2,col = "#E7298A")
abline(v=mean(lam_prob_cmc85_fd$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")


legend(0.01,2,bty="n",legend=c("Past","Present","RCP 4.5","RCP 8.5"),lwd=c(2,2,2,2),lty=1,col=c("#7570B3","#0072B2","#E6AB02","#E7298A"),cex=0.8)

plot(d_past_lambda_2sex, lwd = 2, main = "", xlab = expression("Pr " (lambda[FM]) > 1),,xlim=c(0,1),ylim=c(0,2),col = "#7570B3")
# mtext( "A",side = 3, adj = 0,cex=1.25)
abline(v=mean(geo_prlambda_past$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_lambda_2sex, lwd = 2,col = "#0072B2")
abline(v=mean(geo_prlambda_current$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_cmc45_lambda_2sex, lwd = 2,col = "#E6AB02")
abline(v=mean(geo_prlambdacmc_45$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_cmc85_lambda_2sex, lwd = 2,col = "#E7298A")
abline(v=mean(geo_prlambdacmc_85$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")

dev.off()


pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Densityplot_lambda_ces.pdf",width=4,height=8,useDingbats = F)
par(mar=c(5,5,1.5,0.5),mfrow=c(2,1))
plot(d_past_lambda_fd, lwd = 2, main = "", xlab = expression("Pr " (lambda[F]) > 1),xlim=c(0,1),ylim=c(0,2),col = "#7570B3")
# mtext( "A",side = 3, adj = 0,cex=1.25)
abline(v=mean(lam_prob_past_fd$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_lambda_fd, lwd = 2,col = "#0072B2")
abline(v=mean(lam_prob_current_fd$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_ces45_lambda_fd, lwd = 2,col = "#E6AB02")
abline(v=mean(lam_prob_ces45_fd$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_ces85_lambda_fd, lwd = 2,col = "#E7298A")
abline(v=mean(lam_prob_ces85_fd$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")


legend(0.01,2,bty="n",legend=c("Past","Present","RCP 4.5","RCP 8.5"),lwd=c(2,2,2,2),lty=1,col=c("#7570B3","#0072B2","#E6AB02","#E7298A"),cex=0.8)

plot(d_past_lambda_2sex, lwd = 2, main = "", xlab = expression("Pr " (lambda[FM]) > 1),,xlim=c(0,1),ylim=c(0,2),col = "#7570B3")
# mtext( "A",side = 3, adj = 0,cex=1.25)
abline(v=mean(geo_prlambda_past$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")

lines(d_current_lambda_2sex, lwd = 2,col = "#0072B2")
abline(v=mean(geo_prlambda_current$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")

lines(d_ces45_lambda_2sex, lwd = 2,col = "#E6AB02")
abline(v=mean(geo_prlambdaces_45$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")

lines(d_ces85_lambda_2sex, lwd = 2,col = "#E7298A")
abline(v=mean(geo_prlambdaces_85$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")


dev.off()


pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Densityplot_lambda_GCMs.pdf",width=8,height=7,useDingbats = F)
par(mar=c(5,5,1.5,0.5),mfrow=c(2,2))
plot(d_past_lambda_fd_2sex, lwd = 2, main = "", xlab = expression(paste(Delta,"Pr " (lambda[F-FM]))> 1),xlim=c(-0.3,0.3),ylim=c(0,6),col = "#7570B3")
mtext( "A",side = 3, adj = 0,cex=1.25)
mtext("MIROC 5",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(geo_prlambda_past_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")
lines(d_current_lambda_fd_2sex, lwd = 2,col = "#0072B2")
abline(v=mean(geo_prlambda_current_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")
lines(d_miroc45_lambda_fd_2sex, lwd = 2,col = "#E6AB02")
abline(v=mean(geo_prlambda_miroc45_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")
lines(d_miroc85_lambda_fd_2sex, lwd = 2,col = "#E7298A")
abline(v=mean(geo_prlambda_miroc85_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")

plot(d_past_lambda_fd_2sex, lwd = 2, main = "", xlab = expression(paste(Delta,"Pr " (lambda[F-FM]))> 1),xlim=c(-0.20,0.45),ylim=c(0,6),col = "#7570B3")
mtext( "B",side = 3, adj = 0,cex=1.25)
mtext("ACCESS1-3",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(geo_prlambda_past_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")
lines(d_current_lambda_fd_2sex, lwd = 2,col = "#0072B2")
abline(v=mean(geo_prlambda_current_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")
lines(d_acc45_lambda_fd_2sex, lwd = 2,col = "#E6AB02")
abline(v=mean(geo_prlambda_acc45_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")
lines(d_acc85_lambda_fd_2sex, lwd = 2,col = "#E7298A")
abline(v=mean(geo_prlambda_acc85_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")
legend(0.2,6,bty="n",legend=c("Past","Present","RCP 4.5","RCP 8.5"),lwd=c(2,2,2,2),lty=1,col=c("#7570B3","#0072B2","#E6AB02","#E7298A"),cex=0.8)

plot(d_past_lambda_fd_2sex, lwd = 2, main = "", xlab = expression(paste(Delta,"Pr " (lambda[F-FM]))> 1),xlim=c(-0.20,0.45),ylim=c(0,6),col = "#7570B3")
mtext( "C",side = 3, adj = 0,cex=1.25)
mtext("CMCC-CM",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(geo_prlambda_past_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")
lines(d_current_lambda_fd_2sex, lwd = 2,col = "#0072B2")
abline(v=mean(geo_prlambda_current_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")
lines(d_cmc45_lambda_fd_2sex, lwd = 2,col = "#E6AB02")
abline(v=mean(geo_prlambda_cmc45_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")
lines(d_cmc85_lambda_fd_2sex, lwd = 2,col = "#E7298A")
abline(v=mean(geo_prlambda_cmc85_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")

plot(d_past_lambda_fd_2sex, lwd = 2, main = "", xlab = expression(paste(Delta,"Pr " (lambda[F-FM]))> 1),xlim=c(-0.20,0.45),ylim=c(0,6),col = "#7570B3")
mtext( "D",side = 3, adj = 0,cex=1.25)
mtext("CESM1-BGC",side = 3, adj = 0.5,cex=1,line=0.3)
abline(v=mean(geo_prlambda_past_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#7570B3")
lines(d_current_lambda_fd_2sex, lwd = 2,col = "#0072B2")
abline(v=mean(geo_prlambda_current_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#0072B2")
lines(d_ces45_lambda_fd_2sex, lwd = 2,col = "#E6AB02")
abline(v=mean(geo_prlambda_ces45_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#E6AB02")
lines(d_ces85_lambda_fd_2sex, lwd = 2,col = "#E7298A")
abline(v=mean(geo_prlambda_ces85_fd_twosex$Prlambda,na.rm=TRUE),lty=2,col="#E7298A")

dev.off()
