# Project: Forecasting range shifts of a dioecious plant species under climate change
# Purpose: Build vital rate models (survival, growth,flowering, fertility) for subsequent use in MPMs. 
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

# Demography data ----
poar_allsites <- read.csv("https://www.dropbox.com/s/xk4225mn8btqhbm/demography_allsites.csv?dl=1", stringsAsFactors = F) #common garden data
poar_allsites$census.year<-poar_allsites$year-1 # Add census year to match with climate data
# unique(poar_allsites$census.year)
poar_allsites %>% 
  dplyr::select(everything()) %>% 
  filter(census.year %in% (2015:2016))-> poar_allsites_2015_2016 # Drop the first census to match with the seasonal model (growing and dormant season temperature and precipitation)

# climatic data for the observed period to build the model ----
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

## Observed data (growing and dormant season temperature and precipitation)---- 
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

## Matching seasonal data with observed demographic data ---- 
poar.clim_seasonal <- left_join(x = poar_allsites_2015_2016,y =poar_2015_2016,by=c("site","census.year","Longitude","Latitude")) # merge the demographic data with climatic data for each site
# head(poar.clim_seasonal)

# Mixed effect models (vital rate as a function of size and  climate) -----
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

## Bundle data for model----
data_sites_season <- list( n_sites    = poarclim_seasonal.surv$site %>% n_distinct,
                           n_sources  = poarclim_seasonal.surv$source %>% n_distinct(),
                           # survival data
                           n_blocks_s = poarclim_seasonal.surv$Block %>% n_distinct,
                           site_s   = poarclim_seasonal.surv$site,
                           source_s =  poarclim_seasonal.surv$source,
                           block_s  = poarclim_seasonal.surv$Block,
                           pptgrow_s=as.vector(poarclim_seasonal.surv$z_ppt_t1_grow),
                           pptdorm_s=as.vector(poarclim_seasonal.surv$z_ppt_t1_dorm),
                           tempgrow_s=as.vector(poarclim_seasonal.surv$z_temp_t1_grow),
                           tempdorm_s=as.vector(poarclim_seasonal.surv$z_temp_t1_dorm),
                           site_block_s = data.frame( site_i  = poarclim_seasonal.surv$site,
                                                      block_i = poarclim_seasonal.surv$Block ) %>%
                             unique %>% .$site_i,
                           male_s   = poarclim_seasonal.surv$Sex-1,
                           size_s   = poarclim_seasonal.surv$log_size_t0,
                           y_s      = poarclim_seasonal.surv$surv_t1,
                           n_s      = nrow(poarclim_seasonal.surv),
                           
                           # growth data
                           n_blocks_g = poarclim_seasonal.grow$Block %>% n_distinct,
                           site_g   = poarclim_seasonal.grow$site,
                           source_g =  poarclim_seasonal.grow$source,
                           block_g  = poarclim_seasonal.grow$Block,
                           pptgrow_g=as.vector(poarclim_seasonal.grow$z_ppt_t1_grow),
                           pptdorm_g=as.vector(poarclim_seasonal.grow$z_ppt_t1_dorm),
                           tempgrow_g=as.vector(poarclim_seasonal.grow$z_temp_t1_grow),
                           tempdorm_g=as.vector(poarclim_seasonal.grow$z_temp_t1_dorm),
                           site_block_g = data.frame( site_i  = poarclim_seasonal.grow$site,
                                                      block_i = poarclim_seasonal.grow$Block ) %>%
                             unique %>% .$site_i,
                           male_g   = poarclim_seasonal.grow$Sex-1,
                           size_g   = poarclim_seasonal.grow$log_size_t0,
                           y_g      = poarclim_seasonal.grow$tillerN_t1,
                           n_g      = nrow(poarclim_seasonal.grow),
                           
                           # flowering data
                           n_blocks_f = poarclim_seasonal.flow$Block %>% n_distinct,
                           site_f   = poarclim_seasonal.flow$site,
                           source_f =  poarclim_seasonal.flow$source,
                           block_f  = poarclim_seasonal.flow$Block,
                           pptgrow_f=as.vector(poarclim_seasonal.flow$z_ppt_t1_grow),
                           pptdorm_f=as.vector(poarclim_seasonal.flow$z_ppt_t1_dorm),
                           tempgrow_f=as.vector(poarclim_seasonal.flow$z_temp_t1_grow),
                           tempdorm_f=as.vector(poarclim_seasonal.flow$z_temp_t1_dorm),
                           site_block_f = data.frame( site_i  = poarclim_seasonal.flow$site,
                                                      block_i = poarclim_seasonal.flow$Block ) %>%
                             unique %>% .$site_i,
                           male_f   = poarclim_seasonal.flow$Sex-1,
                           size_f   = poarclim_seasonal.flow$log_size_t1,
                           y_f      = poarclim_seasonal.flow$flow_t1,
                           n_f      = nrow(poarclim_seasonal.flow),
                           
                           # panicle data
                           n_blocks_p = poarclim_seasonal.panic$Block %>% n_distinct,
                           site_p   = poarclim_seasonal.panic$site,
                           source_p =  poarclim_seasonal.panic$source,
                           block_p  = poarclim_seasonal.panic$Block,
                           pptgrow_p=as.vector(poarclim_seasonal.panic$z_ppt_t1_grow),
                           pptdorm_p=as.vector(poarclim_seasonal.panic$z_ppt_t1_dorm),
                           tempgrow_p=as.vector(poarclim_seasonal.panic$z_temp_t1_grow),
                           tempdorm_p=as.vector(poarclim_seasonal.panic$z_temp_t1_dorm),
                           site_block_p = data.frame( site_i  = poarclim_seasonal.panic$site,
                                                      block_i = poarclim_seasonal.panic$Block ) %>%
                             unique %>% .$site_i,
                           male_p   = poarclim_seasonal.panic$Sex-1,
                           size_p   = poarclim_seasonal.panic$log_size_t1,
                           y_p      = poarclim_seasonal.panic$panic_t1,
                           n_p      = nrow(poarclim_seasonal.panic),
                           # viability
                           n_v       = nrow(viab),
                           y_v       = viab$y_viab,
                           tot_seeds_v = viab$tot_seeds_viab,
                           SR_v        = viab$SR,
                           
                           # germination
                           n_m       = nrow(germ),
                           y_m       = germ$y_germ,
                           tot_seeds_m = germ$tot_seeds_germ,
                           SR_m        = germ$SR,
                           
                           # seeds per panicle
                           n_d = nrow(seeds),
                           y_d = seeds$SeedN
                           
)

## Running the stan model----
# sim_pars <- list(
#   warmup = 1000,
#   iter = 4000,
#   thin = 3,
#   chains = 3
# )

# fit the "full" model (poar_season.stan)
# This will take a long time! 
# If you want to run it, uncomment. 
# Otherwise, skip to loading stan output

# fit_allsites_season <- stan(
#  file = "Analysis/stan/poar_season.stan",
#  data = data_sites_season,
#  warmup = sim_pars$warmup,
#  iter = sim_pars$iter,
#  thin = sim_pars$thin,
#  chains = sim_pars$chains)
# saveRDS(fit_allsites_season, 'C:/YOURDIRECTORY')

#load stan output -- this will also take a while, but not as long as running the model from scratch
fit_allsites_season <- readRDS(url("https://www.dropbox.com/scl/fi/1sbhjorfv15ianvdwli96/fit_1.5_0.5.rds?rlkey=rm9sxy29qzc7himps7n6y0xw6&dl=1")) #ignore warning from readRDS

## Chains convergence ----
traceplot(fit_allsites_season,inc_warmup=TRUE, pars = quote_bare(b0_s,bsizesex_s,btempgrowpptgrow_s,bsex_s,bpptdorm_s,
                                                 bsize_s,bpptgrow_s,btempgrow_s,btempdorm_s,
                                                 bpptgrowsex_s,btempgrowsex_s,bpptdormsex_s,
                                                 btempgrowpptgrowsex_s,btempdormpptdormsex_s,bpptgrow2_s,
                                                 btempgrow2_s,bpptgrow2sex_s,bpptdorm2_s,bpptdorm2sex_s,btempdorm2sex_s,
                                                 btempgrow2sex_s))+theme_bw()

traceplot(fit_allsites_season,inc_warmup=TRUE, pars = quote_bare(b0_g,bsizesex_g,btempgrowpptgrow_g,bsex_g,bpptdorm_g,
                                                 bsize_g,bpptgrow_g,btempgrow_g,btempdorm_g,
                                                 bpptgrowsex_g,btempgrowsex_g,bpptdormsex_g,
                                                 btempgrowpptgrowsex_g,btempdormpptdormsex_g,bpptgrow2_g,
                                                 btempgrow2_g,bpptgrow2sex_g,bpptdorm2_g,bpptdorm2sex_g,btempdorm2sex_g,
                                                 btempgrow2sex_g),ncol=4)+theme_bw()

traceplot(fit_allsites_season,inc_warmup=TRUE, pars = quote_bare(b0_f,bsizesex_f,btempgrowpptgrow_f,bsex_f,bpptdorm_f,
                                                 bsize_f,bpptgrow_f,btempgrow_f,btempdorm_f,
                                                 bpptgrowsex_f,btempgrowsex_f,bpptdormsex_f,
                                                 btempgrowpptgrowsex_f,btempdormpptdormsex_f,bpptgrow2_f,
                                                 btempgrow2_f,bpptgrow2sex_f,bpptdorm2_f,bpptdorm2sex_f,btempdorm2sex_f,
                                                 btempgrow2sex_f))+theme_bw()
traceplot(fit_allsites_season, pars = quote_bare(b0_p,bsizesex_p,btempgrowpptgrow_p,bsex_p,bpptdorm_p,
                                                 bsize_p,bpptgrow_p,btempgrow_p,btempdorm_p,
                                                 bpptgrowsex_p,btempgrowsex_p,bpptdormsex_p,
                                                 btempgrowpptgrowsex_p,btempdormpptdormsex_p,bpptgrow2_p,
                                                 btempgrow2_p,bpptgrow2sex_p,bpptdorm2_p,bpptdorm2sex_p,btempdorm2sex_p,
                                                 btempgrow2sex_p))+theme_bw()

## Posterior predictive check for all vital rates ----
#pull out parameter posterior distributions
predS <- rstan::extract(fit_allsites_season, pars = c("predS"))$predS
predG <- rstan::extract(fit_allsites_season, pars = c("predG"))$predG
sigmaG <- rstan::extract(fit_allsites_season, pars = c("sigma"))$sigma
predF <- rstan::extract(fit_allsites_season, pars = c("predF"))$predF
predP <- rstan::extract(fit_allsites_season, pars = c("predP"))$predP
phi_P <- rstan::extract(fit_allsites_season, pars = c("phi_p"))$phi_p
predV <- rstan::extract(fit_allsites_season, pars = c("predV"))$predV
phi_V <- rstan::extract(fit_allsites_season, pars = c("phi_v"))$phi_v
predM <- rstan::extract(fit_allsites_season, pars = c("predM"))$predM
phi_M <- rstan::extract(fit_allsites_season, pars = c("phi_m"))$phi_m
#draw 500 random samples from the joint posterior
n_post_draws <- 500
post_draws <- sample.int(dim(predS)[1], n_post_draws)
#set up simulation output
y_s_sim <- matrix(NA,n_post_draws,length(data_sites_season$y_s))
y_g_sim <- matrix(NA,n_post_draws,length(data_sites_season$y_g))
y_f_sim <- matrix(NA,n_post_draws,length(data_sites_season$y_f))
y_p_sim <- matrix(NA,n_post_draws,length(data_sites_season$y_p))
y_v_sim <- matrix(NA,n_post_draws,length(data_sites_season$y_v))
y_m_sim <- matrix(NA,n_post_draws,length(data_sites_season$y_m))
#loop over the posterior and generate new observations
for(i in 1:n_post_draws){
  print(i)
  ## sample survival data (bernoulli)
  y_s_sim[i,] <- rbinom(n=length(data_sites_season$y_s), size=1, prob = invlogit(predS[i,]))
  ## sample flowering data (bernoulli)
  y_f_sim[i,] <- rbinom(n=length(data_sites_season$y_f), size=1, prob = invlogit(predF[i,]))
  ## sample viability data (beta-binomial)
  y_v_sim[i,] <- rbetabinom(n=length(data_sites_season$y_v), size=data_sites_season$tot_seeds_v, m=predV[i,], s=phi_V[i])
  ## sample germination data (beta-binomial)
  y_m_sim[i,] <- rbetabinom(n=length(data_sites_season$y_m), size=data_sites_season$tot_seeds_m, m=predM[i,], s=phi_M[i])
  ## sample growth data (zero-truncated PIG) -- generates data as legit PIG
  for(j in 1:length(data_sites_season$y_g)){
    ## the pig function appears numerically unstable at low probabilities, so here is a hacky solution
    pig<-dpoisinvgauss(0:1000,mean=predG[i,j],shape=(sigmaG[i]*predG[i,j]))
    pig[is.nan(pig) | is.infinite(pig)] <- 0
    pig_trunc_prob <- pig[2:1001] / (1 - ((1 - sum(pig)) + pig[1]))
    y_g_sim[i,j] <- sample(x=1:1000,size=1,replace=T,prob=pig_trunc_prob)
  } 
  ## sample panicle data (zero-truncated NB)
  for(j in 1:length(data_sites_season$y_p)){
    y_p_sim[i,j] <- sample(x=1:1000,size=1,replace=T,prob=dnbinom(1:1000, mu = exp(predP[i,j]), size=phi_P[i]) / (1 - dnbinom(0, mu = exp(predP[i,j]), size=phi_P[i])))
  }
}

# plot ppc overlay
ppc_dens_overlay(data_sites_season$y_s, na.omit(y_s_sim))+
  xlab("Survival status")+ylab("Density")+
  ggtitle(("Survival"))+theme(legend.position = "none")->ppc_surv
ppc_dens_overlay(data_sites_season$y_g, na.omit(y_g_sim))+xlim(0, 50)+
  xlab("Number of tillers")+ylab("Density")+
  ggtitle(("Growth"))+theme(legend.position = "none")->ppc_grow
ppc_dens_overlay(data_sites_season$y_f, y_f_sim)+
  xlab("Flowering status")+ylab("Density")+
  ggtitle(("Flowering"))+theme(legend.position = "none")->ppc_flow
ppc_dens_overlay(data_sites_season$y_p, y_p_sim)+xlim(0, 30)+
  xlab("Number of panicles")+ylab("Density")+
  ggtitle(("Panicles"))+theme(legend.position = "none")->ppc_panic
ppc_dens_overlay(data_sites_season$y_v, y_v_sim)+
  xlab("Number of viable seeds")+ylab("Density")+
  ggtitle(("Seed viability"))+theme(legend.position = "none")->ppc_viab
ppc_dens_overlay(data_sites_season$y_m, y_m_sim)+
  xlab("Number of germinants")+ylab("Density")+
  ggtitle(("Seed germination"))+theme(legend.position = "none")->ppc_germ

#print figure
# pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/PPC.pdf",useDingbats = F,height=10,width=6)
# multiplot(ppc_surv,ppc_flow,ppc_viab,
#           ppc_grow,ppc_panic,ppc_germ,cols=2)
# dev.off()

## Posterior mean values for each vital rate -----
posterior_season <- as.array(fit_allsites_season)
color_scheme_set("orange")
surv_clim<-mcmc_intervals(posterior_season, pars = quote_bare(b0_s,bsizesex_s,btempgrowpptgrow_s,btempdormpptdorm_s,bsex_s,bpptdorm_s,
                                                              bsize_s,bpptgrow_s,btempgrow_s,btempdorm_s,
                                                              bpptgrowsex_s,btempgrowsex_s,bpptdormsex_s,btempdormsex_s,
                                                              btempgrowpptgrowsex_s,btempdormpptdormsex_s,bpptgrow2_s,
                                                              btempgrow2_s,btempdorm2_s,bpptgrow2sex_s,bpptdorm2_s,bpptdorm2sex_s,btempdorm2sex_s,
                                                              btempgrow2sex_s)) + 
  ggplot2::scale_y_discrete(limits = c("b0_s","bsex_s", "bsize_s","bsizesex_s",
                                       "bpptgrow_s","bpptdorm_s","btempgrow_s","btempdorm_s",
                                       "bpptgrowsex_s","bpptdormsex_s", "btempdormsex_s","btempgrowsex_s",
                                       "btempgrowpptgrowsex_s", "btempdormpptdormsex_s", 
                                       "btempgrowpptgrow_s","btempdormpptdorm_s",
                                       "bpptgrow2_s","bpptdorm2_s","btempdorm2_s","btempgrow2_s",
                                       "bpptdorm2sex_s","bpptgrow2sex_s","btempdorm2sex_s", "btempgrow2sex_s"),
                            labels=c("b0_s"="Intercept",
                                     "bsex_s"="sex",
                                     "bsize_s"="size",
                                     "bsizesex_s"="size:sex",
                                     "bpptgrow_s"="pptgrow",
                                     "bpptdorm_s"="pptdorm",
                                     "btempgrow_s"="tempgrow",
                                     "btempdorm_s"="tempdorm",
                                     "bpptgrowsex_s"="pptgrow:sex",
                                     "bpptdormsex_s"="pptdorm:sex",
                                     "btempdormsex_s"="tempdorm:sex",
                                     "btempgrowsex_s"="tempgrow:sex",  
                                     "btempgrowpptgrowsex_s"="tempgrow:pptgrow:sex",
                                     "btempdormpptdormsex_s"="tempdorm:pptdorm:sex",
                                     "btempgrowpptgrow_s"="tempgrow:pptgrow",
                                     "btempdormpptdorm_s"="tempdorm:pptdorm",
                                     "bpptgrow2_s"=expression(pptgrow^2),
                                     "btempgrow2_s"=expression(tempgrow^2),
                                     "bpptdorm2_s"=expression(pptdorm^2),
                                     "btempdorm2_s"=expression(tempdorm^2),
                                     "bpptdorm2sex_s"=expression(pptdorm^2:sex),
                                     "btempdorm2sex_s"=expression(tempdorm^2:sex),
                                     "bpptgrow2sex_s"=expression(pptgrow^2:sex),
                                     "btempgrow2sex_s"=expression(tempgrow^2:sex)))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Survival)")+
  xlim(-2,2)+
  ggtitle("A") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

grow_clim<-mcmc_intervals(posterior_season, pars = quote_bare(b0_g,bsizesex_g,btempgrowpptgrow_g,btempdormpptdorm_g,bsex_g,bpptdorm_g,
                                                              bsize_g,bpptgrow_g,btempgrow_g,btempdorm_g,
                                                              bpptgrowsex_g,btempgrowsex_g,bpptdormsex_g,btempdormsex_g,
                                                              btempgrowpptgrowsex_g,btempdormpptdormsex_g,bpptgrow2_g,
                                                              btempgrow2_g,btempdorm2_g,bpptgrow2sex_g,bpptdorm2_g,bpptdorm2sex_g,btempdorm2sex_g,
                                                              btempgrow2sex_g)) + 
  ggplot2::scale_y_discrete(limits = c("b0_g","bsex_g", "bsize_g","bsizesex_g",
                                       "bpptgrow_g","bpptdorm_g","btempgrow_g","btempdorm_g",
                                       "bpptgrowsex_g","bpptdormsex_g", "btempdormsex_g","btempgrowsex_g",
                                       "btempgrowpptgrowsex_g", "btempdormpptdormsex_g", 
                                       "btempgrowpptgrow_g","btempdormpptdorm_g",
                                       "bpptgrow2_g","bpptdorm2_g","btempdorm2_g","btempgrow2_g",
                                       "bpptdorm2sex_g","bpptgrow2sex_g","btempdorm2sex_g", "btempgrow2sex_g"),
                            labels=c("b0_g"="Intercept",
                                     "bsize_g"="size",
                                     "bsex_g"="sex",
                                     "bsizesex_g"="size:sex",
                                     "bpptgrow_g"="pptgrow",
                                     "bpptdorm_g"="pptdorm",
                                     "btempgrow_g"="tempgrow",
                                     "btempdorm_g"="tempdorm",
                                     "bpptgrowsex_g"="pptgrow:sex",
                                     "bpptdormsex_g"="pptdorm:sex",
                                     "btempdormsex_g"="tempdorm:sex",
                                     "btempgrowsex_g"="tempgrow:sex",  
                                     "btempgrowpptgrowsex_g"="tempgrow:pptgrow:sex",
                                     "btempdormpptdormsex_g"="tempdorm:pptdorm:sex",
                                     "btempgrowpptgrow_g"="tempgrow:pptgrow",
                                     "btempdormpptdorm_g"="tempdorm:pptdorm",
                                     "bpptgrow2_g"=expression(pptgrow^2),
                                     "btempgrow2_g"=expression(tempgrow^2),
                                     "bpptdorm2_g"=expression(pptdorm^2),
                                     "btempdorm2_g"=expression(tempdorm^2),
                                     "bpptdorm2sex_g"=expression(pptdorm^2:sex),
                                     "btempdorm2sex_g"=expression(tempdorm^2:sex),
                                     "bpptgrow2sex_g"=expression(pptgrow^2:sex),
                                     "btempgrow2sex_g"=expression(tempgrow^2:sex)))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Growth)")+
  ylab(NULL) +
  xlim(-1.5,1.5)+
  ggtitle("B") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

flow_clim<-mcmc_intervals(posterior_season, pars = quote_bare(b0_g,bsizesex_f,btempgrowpptgrow_f,btempdormpptdorm_f,bsex_f,bpptdorm_f,
                                                              bsize_f,bpptgrow_f,btempgrow_f,btempdorm_f,
                                                              bpptgrowsex_f,btempgrowsex_f,bpptdormsex_f,btempdormsex_f,
                                                              btempgrowpptgrowsex_f,btempdormpptdormsex_f,bpptgrow2_f,
                                                              btempgrow2_f,btempdorm2_f,bpptgrow2sex_f,bpptdorm2_f,bpptdorm2sex_f,btempdorm2sex_f,
                                                              btempgrow2sex_f)) + 
  ggplot2::scale_y_discrete(limits = c("b0_f","bsex_f", "bsize_f","bsizesex_f",
                                       "bpptgrow_f","bpptdorm_f","btempgrow_f","btempdorm_f",
                                       "bpptgrowsex_f","bpptdormsex_f", "btempdormsex_f","btempgrowsex_f",
                                       "btempgrowpptgrowsex_f", "btempdormpptdormsex_f", 
                                       "btempgrowpptgrow_f","btempdormpptdorm_f",
                                       "bpptgrow2_f","bpptdorm2_f","btempdorm2_f","btempgrow2_f",
                                       "bpptdorm2sex_f","bpptgrow2sex_f","btempdorm2sex_f", "btempgrow2sex_f"),
                            labels=c("b0_f"="Intercept",
                                     "bsize_f"="size",
                                     "bsex_f"="sex",
                                     "bsizesex_f"="size:sex",
                                     "bpptgrow_f"="pptgrow",
                                     "bpptdorm_f"="pptdorm",
                                     "btempgrow_f"="tempgrow",
                                     "btempdorm_f"="tempdorm",
                                     "bpptgrowsex_f"="pptgrow:sex",
                                     "bpptdormsex_f"="pptdorm:sex",
                                     "btempdormsex_f"="tempdorm:sex",
                                     "btempgrowsex_f"="tempgrow:sex",  
                                     "btempgrowpptgrowsex_f"="tempgrow:pptgrow:sex",
                                     "btempdormpptdormsex_f"="tempdorm:pptdorm:sex",
                                     "btempgrowpptgrow_f"="tempgrow:pptgrow",
                                     "btempdormpptdorm_f"="tempdorm:pptdorm",
                                     "bpptgrow2_f"=expression(pptgrow^2),
                                     "btempgrow2_f"=expression(tempgrow^2),
                                     "bpptdorm2_f"=expression(pptdorm^2),
                                     "btempdorm2_f"=expression(tempdorm^2),
                                     "bpptdorm2sex_f"=expression(pptdorm^2:sex),
                                     "btempdorm2sex_f"=expression(tempdorm^2:sex),
                                     "bpptgrow2sex_f"=expression(pptgrow^2:sex),
                                     "btempgrow2sex_f"=expression(tempgrow^2:sex)))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "black") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Flowering)")+
  ylab(NULL) +
  xlim(-2,2)+
  ggtitle("C") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

panic_clim<-mcmc_intervals(posterior_season, pars = quote_bare(b0_p,bsizesex_p,btempgrowpptgrow_p,btempdormpptdorm_p,bsex_p,bpptdorm_p,
                                                               bsize_p,bpptgrow_p,btempgrow_p,btempdorm_p,
                                                               bpptgrowsex_p,btempgrowsex_p,bpptdormsex_p,btempdormsex_p,
                                                               btempgrowpptgrowsex_p,btempdormpptdormsex_p,bpptgrow2_p,
                                                               btempgrow2_p,btempdorm2_p,bpptgrow2sex_p,bpptdorm2_p,bpptdorm2sex_p,btempdorm2sex_p,
                                                               btempgrow2sex_p)) + 
  ggplot2::scale_y_discrete(limits = c("b0_p","bsex_p", "bsize_p","bsizesex_p",
                                       "bpptgrow_p","bpptdorm_p","btempgrow_p","btempdorm_p",
                                       "bpptgrowsex_p","bpptdormsex_p", "btempdormsex_p","btempgrowsex_p",
                                       "btempgrowpptgrowsex_p", "btempdormpptdormsex_p", 
                                       "btempgrowpptgrow_p","btempdormpptdorm_p",
                                       "bpptgrow2_p","bpptdorm2_p","btempdorm2_p","btempgrow2_p",
                                       "bpptdorm2sex_p","bpptgrow2sex_p","btempdorm2sex_p", "btempgrow2sex_p"),
                            labels=c("b0_p"="Intercept","bsize_p"="size","bsex_p"="sex",
                                     "bsizesex_p"="size:sex",
                                     "bpptgrow_p"="pptgrow","bpptdorm_p"="pptdorm","btempgrow_p"="tempgrow","btempdorm_p"="tempdorm",
                                     "bpptgrowsex_p"="pptgrow:sex","bpptdormsex_p"="pptdorm:sex","btempdormsex_p"="tempdorm:sex","btempgrowsex_p"="tempgrow:sex",  
                                     "btempgrowpptgrowsex_p"="tempgrow:pptgrow:sex","btempdormpptdormsex_p"="tempdorm:pptdorm:sex",
                                     "btempgrowpptgrow_p"="tempgrow:pptgrow", "btempdormpptdorm_p"="tempdorm:pptdorm",
                                     "bpptgrow2_p"=expression(pptgrow^2),"btempgrow2_p"=expression(tempgrow^2),"bpptdorm2_p"=expression(pptdorm^2),"btempdorm2_p"=expression(tempdorm^2),
                                     "bpptdorm2sex_p"=expression(pptdorm^2:sex),"btempdorm2sex_p"=expression(tempdorm^2:sex),"bpptgrow2sex_p"=expression(pptgrow^2:sex),"btempgrow2sex_p"=expression(tempgrow^2:sex)))+
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, alpha = 0.6, color = "white") +
  labs(color = "Interaction type:")+
  xlab("Posterior estimates (Panicles)")+
  ylab(NULL) +
  xlim(-1,1)+
  ggtitle("D") +
  # geom_rect(xmin = 0, xmax=2.25, ymin = 0, ymax = 25, alpha = 0.006, fill = "#F4B400", color = NA)+
  # geom_rect(xmin = -2.25, xmax = 0, ymin = 0, ymax = 25, alpha = 0.006, fill = "#0F9D58", color = NA)+
  theme_pubr()+
  theme(axis.title.x = element_text(family = "Helvetica",colour="black", size = 13),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y  = element_text(size = 10),
        axis.line.x = element_line(linewidth = 0.1),
        axis.line.y = element_line(linewidth = 0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))

# pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Posterior_mean.pdf",useDingbats = F,height=12,width=10)
# ggarrange(surv_clim,grow_clim,flow_clim,panic_clim + rremove("ylab"), ncol = 2, nrow = 2)
# dev.off()

# Plot vital rate figure (2D) ----
poar <- poar.clim_seasonal  
fit_full <- fit_allsites_season  

# bin size groups for visualization
size_bin_num <- 1

## re-format data for plotting
## Survival
poar %>% 
  subset( tillerN_t0 > 0 )%>%
  dplyr::select(census.year, Code, site, Block, Sex, 
                Longitude, Latitude, 
                tillerN_t0, surv_t1,zpptgrow,zpptdorm,ztempgrow,ztempdorm,site)%>% 
  na.omit %>% 
  mutate( site         = site %>% as.factor %>% as.numeric,
          Block = Block %>% as.factor %>% as.numeric,
          Sex          = Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric ) %>%
  mutate( log_size_t0   = log(tillerN_t0),
          z_size_t0   = tillerN_t0 %>% scale %>% .[,1],
          z_ppt_t1_grow = zpptgrow,
          z_ppt_t1_dorm = zpptdorm ,
          z_temp_t1_grow = ztempgrow ,
          z_temp_t1_dorm = ztempdorm)->poar.surv


poar_surv_binned <- poar.surv %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t0,size_bin_num))) %>% 
  group_by(site,census.year,Sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t0),
            mean_surv = mean(surv_t1),
            pptgrow = unique(z_ppt_t1_grow),
            pptdorm = unique(z_ppt_t1_dorm),
            tempgrow = unique(z_temp_t1_grow),
            tempdorm = unique(z_temp_t1_dorm),
            bin_n = n())

surv_mean_sizes <- poar_surv_binned %>% group_by(Sex,size_bin) %>% summarise(size = mean(mean_size))
surv_mean_pptgrow <- poar_surv_binned %>% group_by(Sex,size_bin) %>% summarise(pptgrow = mean(pptgrow))
surv_mean_tempgrow <- poar_surv_binned %>% group_by(Sex,size_bin) %>% summarise(tempgrow = mean(tempgrow))
surv_mean_pptdorm <- poar_surv_binned %>% group_by(Sex,size_bin) %>% summarise(pptdorm = mean(pptdorm))
surv_mean_tempdorm <- poar_surv_binned %>% group_by(Sex,size_bin) %>% summarise(tempdorm = mean(tempdorm))

## Growth
poar %>% 
  subset( tillerN_t0 > 0 & tillerN_t1 > 0)%>%
  dplyr::select( census.year, Code, site, Block, Sex, 
                 Longitude, Latitude, 
                 tillerN_t0, tillerN_t1,surv_t1,zpptgrow,zpptdorm,ztempgrow,ztempdorm,site)%>% 
  na.omit %>% 
  mutate( site         = site %>% as.factor %>% as.numeric,
          Block = Block %>% as.factor %>% as.numeric,
          Sex          = Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric ) %>%
  mutate( log_size_t0   = log(tillerN_t0),
          z_size_t0   = tillerN_t0,
          z_ppt_t1_grow = zpptgrow,
          z_ppt_t1_dorm = zpptdorm,
          z_temp_t1_grow = ztempgrow,
          z_temp_t1_dorm = ztempdorm)->poar.grow

poar_grow_binned <- poar.grow %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t0,size_bin_num))) %>% 
  group_by(site,census.year,Sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t0),
            mean_grow = mean(tillerN_t1),
            pptgrow = unique(z_ppt_t1_grow),
            pptdorm = unique(z_ppt_t1_dorm),
            tempgrow = unique(z_temp_t1_grow),
            tempdorm = unique(z_temp_t1_dorm),
            bin_n = n())
grow_mean_sizes <- poar_grow_binned %>% group_by(Sex,size_bin) %>% summarise(size = mean(mean_size))
grow_mean_pptgrow <- poar_grow_binned %>% group_by(Sex,size_bin) %>% summarise(pptgrow = mean(pptgrow))
grow_mean_tempgrow <- poar_grow_binned %>% group_by(Sex,size_bin) %>% summarise(tempgrow = mean(tempgrow))
grow_mean_pptdorm <- poar_grow_binned %>% group_by(Sex,size_bin) %>% summarise(pptdorm = mean(pptdorm))
grow_mean_tempdorm <- poar_grow_binned %>% group_by(Sex,size_bin) %>% summarise(tempdorm = mean(tempdorm))

## Flowering
poar %>% 
  subset( tillerN_t1 > 0 )%>%
  dplyr::select(census.year, Code, site, Block, Sex, 
                Longitude, Latitude, 
                tillerN_t1, flowerN_t1,flow_t1,zpptgrow,zpptdorm,ztempgrow,ztempdorm,site)%>% 
  na.omit %>% 
  mutate( site         = site %>% as.factor %>% as.numeric,
          Block = Block %>% as.factor %>% as.numeric,
          Sex          = Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric ) %>%
  mutate( log_size_t1   = log(tillerN_t1),
          z_size_t1 = tillerN_t1 %>% scale %>% .[,1],
          z_ppt_t1_grow = zpptgrow,
          z_ppt_t1_dorm = zpptdorm ,
          z_temp_t1_grow = ztempgrow,
          z_temp_t1_dorm = ztempdorm)->poar.flow

poar_flow_binned <- poar.flow %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t1,size_bin_num))) %>% 
  group_by(site,census.year,Sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t1),
            mean_flow = mean(flow_t1),
            pptgrow = unique(z_ppt_t1_grow),
            pptdorm = unique(z_ppt_t1_dorm),
            tempgrow = unique(z_temp_t1_grow),
            tempdorm = unique(z_temp_t1_dorm),
            bin_n = n())
flow_mean_sizes <- poar_flow_binned %>% group_by(Sex,size_bin) %>% summarise(size = mean(mean_size))
flow_mean_pptgrow <- poar_flow_binned %>% group_by(Sex,size_bin) %>% summarise(pptgrow = mean(pptgrow))
flow_mean_tempgrow <- poar_flow_binned %>% group_by(Sex,size_bin) %>% summarise(tempgrow = mean(tempgrow))
flow_mean_pptdorm <- poar_flow_binned %>% group_by(Sex,size_bin) %>% summarise(pptdorm = mean(pptdorm))
flow_mean_tempdorm <- poar_flow_binned %>% group_by(Sex,size_bin) %>% summarise(tempdorm = mean(tempdorm))

## Panicles
poar %>% 
  subset( flowerN_t1 > 0 & tillerN_t1 > 0 )%>%
  dplyr::select( census.year, Code, site, Block, Sex, 
                 Longitude, Latitude, 
                 tillerN_t1, flowerN_t1,zpptgrow,zpptdorm,ztempgrow,ztempdorm,site)%>% 
  na.omit %>% 
  mutate( site         = site %>% as.factor %>% as.numeric,
          Block = Block %>% as.factor %>% as.numeric,
          Sex          = Sex %>% as.factor %>% as.numeric,
          source = Code %>% as.factor %>% as.numeric,
          panic_t1 = flowerN_t1) %>%
  mutate( 
    log_size_t1   = log(tillerN_t1),
    z_size_t1 = tillerN_t1 %>% scale %>% .[,1],
    z_ppt_t1_grow = zpptgrow ,
    z_ppt_t1_dorm = zpptdorm ,
    z_temp_t1_grow = ztempgrow ,
    z_temp_t1_dorm = ztempdorm)->poar.panic

poar_panic_binned <- poar.panic %>% 
  mutate(size_bin = as.integer(cut_number(log_size_t1,size_bin_num))) %>% 
  group_by(site,census.year,Sex,size_bin) %>% 
  summarise(mean_size = mean(log_size_t1),
            mean_panic = mean(panic_t1),
            pptgrow = unique(z_ppt_t1_grow),
            pptdorm = unique(z_ppt_t1_dorm),
            tempgrow = unique(z_temp_t1_grow),
            tempdorm = unique(z_temp_t1_dorm),
            bin_n = n())

panic_mean_sizes <- poar_panic_binned %>% group_by(Sex,size_bin) %>% summarise(size = mean(mean_size))
panic_mean_pptgrow <- poar_panic_binned %>% group_by(Sex,size_bin) %>% summarise(pptgrow = mean(pptgrow))
panic_mean_tempgrow <- poar_panic_binned %>% group_by(Sex,size_bin) %>% summarise(tempgrow = mean(tempgrow))
panic_mean_pptdorm <- poar_panic_binned %>% group_by(Sex,size_bin) %>% summarise(pptdorm = mean(pptdorm))
panic_mean_tempdorm <- poar_panic_binned %>% group_by(Sex,size_bin) %>% summarise(tempdorm = mean(tempdorm))

# pull out stan coefficients
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

## Precipitation of the growing season 
# sample vital rate functions from the join posterior
pptgrow_seq <- seq(min(poar_surv_binned$pptgrow),max(poar_surv_binned$pptgrow),length.out=30)
n_post_draws <- 2000
## Temperature of the growing season 
# sample vital rate functions from the join posterior
tempgrow_seq <- seq(min(poar_surv_binned$tempgrow),max(poar_surv_binned$tempgrow),length.out=30)
## Precipitation of the dormant season
# sample vital rate functions from the join posterior
pptdorm_seq <- seq(min(poar_surv_binned$pptdorm),max(poar_surv_binned$pptdorm),length.out=30)
## Temperature of the dormant season
# sample vital rate functions from the join posterior
tempdorm_seq <- seq(min(poar_surv_binned$tempdorm),max(poar_surv_binned$tempdorm),length.out=30)

# set up figure
sex_cols <- RColorBrewer::brewer.pal(8, "Dark2")[c(2,1)]
sex_lty <- c(1,1)
bin_shapes <- 15
diff_col <- RColorBrewer::brewer.pal(8, "Dark2")[6]
diff_alpha <- 0.35
graph<-list("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")
layout.matrix <- rbind(matrix(1:4, nrow = 1, ncol = 4, byrow = F),
                       matrix(5:8, nrow = 1, ncol = 4, byrow = F),
                       matrix(9:12, nrow = 1, ncol = 4, byrow = F),
                       matrix(13:16, nrow = 1, ncol = 4, byrow = F))
#print figure
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/vital_rates_v1.pdf",height = 12,width = 14,useDingbats = F)
layout(mat = layout.matrix,
       heights = c(2, 2, 2,2), # Heights of one row
       widths = c(2, 2, 2,2))
# layout.show(16)

par(oma=c(4,1,0.5,0.5))
with(poar_surv_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(pptgrow[size_bin==i]*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),mean_surv[size_bin==i],type="n",xaxt="n",ylim=c(0,1),
         xlab=" ",ylab=" ");box()
    if(i==1){mtext("Pr(survival)",side=2,line=3,cex=1.5)}
    mylabel3 <- paste(graph[[1]])
    mtext(mylabel3,side = 3, adj = 0,cex=1.5)
    for(s in 1:2){
      points(pptgrow[Sex==s & size_bin==i]*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),mean_surv[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(pptgrow_seq*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),invlogit(mean(surv_coef$b0_s) + 
                         mean(surv_coef$bsize_s) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i]+
                         mean(surv_coef$bsizesex_s) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                         mean(surv_coef$bsex_s) * (s-1) +
                         mean(surv_coef$bpptgrow_s) * pptgrow_seq +
                         mean(surv_coef$bpptdorm_s) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s]  +
                         mean(surv_coef$btempgrow_s) * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] +
                         mean(surv_coef$btempdorm_s) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] +
                         mean(surv_coef$bpptgrowsex_s) * pptgrow_seq * (s-1) +
                         mean(surv_coef$bpptdormsex_s) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s]  * (s-1) +
                         mean(surv_coef$btempgrowsex_s) * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] * (s-1) +
                        mean(surv_coef$btempdormsex_s) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * (s-1) +
                         mean(surv_coef$btempdormpptdorm_s) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] +
                         mean(surv_coef$btempgrowpptgrow_s) * pptgrow_seq * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] +
                         mean(surv_coef$btempdormpptdormsex_s) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] * (s-1) +
                         mean(surv_coef$btempgrowpptgrowsex_s) * pptgrow_seq * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s]  * (s-1) +
                         mean(surv_coef$bpptgrow2_s) * (pptgrow_seq^2) +
                         mean(surv_coef$bpptdorm2_s) * (surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s])^2 +
                         mean(surv_coef$btempgrow2_s) * (surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s])^2 +
                         mean(surv_coef$btempdorm2_s) * (surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s])^2 +
                         mean(surv_coef$bpptgrow2sex_s) * (pptgrow_seq^2) * (s-1) +
                         mean(surv_coef$bpptdorm2sex_s) * (surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s])^2 * (s-1) +
                         mean(surv_coef$btempgrow2sex_s) * (surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s])^2 * (s-1) +
                         mean(surv_coef$btempdorm2sex_s) * (surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s])^2 * (s-1)
    ),lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_surv_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(tempgrow[size_bin==i]*sd(poar_2015_2016$tempgrow)+ mean(poar_2015_2016$tempgrow),mean_surv[size_bin==i],type="n",ylim=c(0,1),
         xlab=" ",ylab=" ",xaxt="n");box()
     mylabel3 <- paste(graph[[2]])
     mtext(mylabel3,side = 3, adj = 0,cex=1.5)
    for(s in 1:2){
      points(tempgrow[Sex==s & size_bin==i]*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow),mean_surv[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(tempgrow_seq*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),invlogit(mean(surv_coef$b0_s) +
                        mean(surv_coef$bsize_s) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i] +
                        mean(surv_coef$bsizesex_s) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                        mean(surv_coef$bsex_s) * (s-1) +
                        mean(surv_coef$bpptgrow_s) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s]  +
                        mean(surv_coef$bpptdorm_s) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s]  +
                        mean(surv_coef$btempgrow_s) * tempgrow_seq +
                        mean(surv_coef$btempdorm_s) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] +
                        mean(surv_coef$bpptgrowsex_s) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * (s-1) +
                        mean(surv_coef$bpptdormsex_s) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s]  * (s-1) +
                        mean(surv_coef$btempdormpptdorm_s) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] +
                        mean(surv_coef$btempgrowpptgrow_s) * tempgrow_seq * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] +
                        mean(surv_coef$btempgrowsex_s) * tempgrow_seq * (s-1) +
                        mean(surv_coef$btempdormsex_s) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * (s-1) +
                        mean(surv_coef$btempdormpptdormsex_s) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] * (s-1) +
                        mean(surv_coef$btempgrowpptgrowsex_s) * tempgrow_seq * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s]  * (s-1) +
                        mean(surv_coef$bpptgrow2_s) * (surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s])^2 +
                        mean(surv_coef$bpptdorm2_s) * (surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s])^2 +
                        mean(surv_coef$btempgrow2_s) * (tempgrow_seq^2) +
                        mean(surv_coef$btempdorm2_s) * (surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s])^2 +
                        mean(surv_coef$bpptgrow2sex_s) * (surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s])^2 * (s-1) +
                        mean(surv_coef$bpptdorm2sex_s) * (surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s])^2 * (s-1) +
                        mean(surv_coef$btempgrow2sex_s) * (tempgrow_seq^2) * (s-1) +
                        mean(surv_coef$btempdorm2sex_s) * (surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s])^2 * (s-1)
    
    ),lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_surv_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(pptdorm[size_bin==i]*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),mean_surv[size_bin==i],type="n",ylim=c(0,1),
         xlab=" ",ylab=" ",xaxt="n");box()
     mylabel3 <- paste(graph[[3]])
     mtext(mylabel3,side = 3, adj = 0,cex=1.2)
    for(s in 1:2){
      points(pptdorm[Sex==s & size_bin==i]*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),mean_surv[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(pptdorm_seq*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),invlogit(mean(surv_coef$b0_s) + 
                         mean(surv_coef$bsize_s) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i]+
                         mean(surv_coef$bsizesex_s) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                         mean(surv_coef$bsex_s) * (s-1) +
                         mean(surv_coef$bpptgrow_s) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] +
                         mean(surv_coef$bpptdorm_s) *  pptdorm_seq +
                         mean(surv_coef$btempgrow_s) * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] +
                         mean(surv_coef$btempdorm_s) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] +
                         mean(surv_coef$bpptgrowsex_s) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * (s-1) +
                         mean(surv_coef$bpptdormsex_s) * pptdorm_seq * (s-1) +
                         mean(surv_coef$btempgrowsex_s) * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] * (s-1) +
                        mean(surv_coef$btempdormsex_s) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * (s-1) +
                         mean(surv_coef$btempdormpptdorm_s) * pptdorm_seq * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] +
                         mean(surv_coef$btempgrowpptgrow_s) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] +
                         mean(surv_coef$btempdormpptdormsex_s) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * pptdorm_seq * (s-1) +
                         mean(surv_coef$btempgrowpptgrowsex_s) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] * (s-1) +
                         mean(surv_coef$bpptgrow2_s) * (surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s])^2 +
                         mean(surv_coef$bpptdorm2_s) * (pptdorm_seq^2) +
                         mean(surv_coef$btempgrow2_s) * (surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s])^2 +
                         mean(surv_coef$btempdorm2_s) * (surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s])^2 +
                         mean(surv_coef$bpptgrow2sex_s) * (surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s])^2 * (s-1) +
                         mean(surv_coef$bpptdorm2sex_s) * (pptdorm_seq^2) * (s-1) +
                         mean(surv_coef$btempgrow2sex_s) * (surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s])^2 * (s-1) +
                         mean(surv_coef$btempdorm2sex_s) * (surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s])^2 * (s-1)
    ),lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
    par(mar=c(2,4,0.5,0))
  }
})

with(poar_surv_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(tempdorm[size_bin==i]*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mean_surv[size_bin==i],type="n",ylim=c(0,1),
         xlab=" ",ylab=" ",xaxt="n");box()
     mylabel3 <- paste(graph[[4]])
     mtext(mylabel3,side = 3, adj = 0,cex=1.5)
    for(s in 1:2){
      points(tempdorm[Sex==s & size_bin==i]*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mean_surv[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),invlogit(mean(surv_coef$b0_s) + 
                         mean(surv_coef$bsize_s) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i] +
                         mean(surv_coef$bsizesex_s) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                         mean(surv_coef$bsex_s) * (s-1) +
                         mean(surv_coef$bpptgrow_s) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] +
                         mean(surv_coef$bpptdorm_s) *  surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] +
                         mean(surv_coef$btempgrow_s) * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] +
                         mean(surv_coef$btempdorm_s) * tempdorm_seq +
                         mean(surv_coef$bpptgrowsex_s) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * (s-1) +
                         mean(surv_coef$bpptdormsex_s) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] * (s-1) +
                         mean(surv_coef$btempgrowsex_s) * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] * (s-1) +
                        mean(surv_coef$btempdormsex_s) * tempdorm_seq * (s-1) +
                         mean(surv_coef$btempdormpptdorm_s) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] +
                         mean(surv_coef$btempgrowpptgrow_s) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] +
                         mean(surv_coef$btempdormpptdormsex_s) * tempdorm_seq * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] * (s-1) +
                         mean(surv_coef$btempgrowpptgrowsex_s) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] * (s-1) +
                         mean(surv_coef$bpptgrow2_s) * (surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s])^2 +
                         mean(surv_coef$bpptdorm2_s) * (surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s])^2 +
                         mean(surv_coef$btempgrow2_s) * (surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s])^2 +
                         mean(surv_coef$btempdorm2_s) * (tempdorm_seq^2) +
                         mean(surv_coef$bpptgrow2sex_s) * (surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s])^2 * (s-1) +
                        mean(surv_coef$bpptdorm2sex_s) * (surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s])^2 * (s-1) +
                        mean(surv_coef$btempgrow2sex_s) * (surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s])^2 * (s-1) +
                        mean(surv_coef$btempdorm2sex_s) * (tempdorm_seq^2) * (s-1)
    ),lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})


with(poar_grow_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(pptgrow[size_bin==i]*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),mean_grow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n",ylim=c(0,50));box()    
    if(i==1){mtext("#tillers",side=2,line=3,cex=1.5)}
    mylabel3 <- paste(graph[[5]])
    mtext(mylabel3,side = 3, adj = 0,cex=1.5)
    for(s in 1:2){
      points(pptgrow[Sex==s & size_bin==i]*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow) ,mean_grow[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2,ylim=c(0,50))
      lines(pptgrow_seq*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow) ,
            exp(mean(grow_coef$b0_g) + 
                         mean(grow_coef$bsize_g) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i]+
                         mean(grow_coef$bsizesex_g) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                         mean(grow_coef$bsex_g) * (s-1) +
                         mean(grow_coef$bpptgrow_g) * pptgrow_seq +
                         mean(grow_coef$bpptdorm_g) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s]  +
                         mean(grow_coef$btempgrow_g) * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] +
                         mean(grow_coef$btempdorm_g) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] +
                         mean(grow_coef$bpptgrowsex_g) * pptgrow_seq * (s-1) +
                         mean(grow_coef$bpptdormsex_g) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s]  * (s-1) +
                         mean(grow_coef$btempgrowsex_g) * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] * (s-1) +
                        mean(grow_coef$btempdormsex_g) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * (s-1) +
                         mean(grow_coef$btempdormpptdorm_g) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] +
                         mean(grow_coef$btempgrowpptgrow_g) * pptgrow_seq * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] +
                         mean(grow_coef$btempdormpptdormsex_g) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] * (s-1) +
                         mean(grow_coef$btempgrowpptgrowsex_g) * pptgrow_seq * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s]  * (s-1) +
                         mean(grow_coef$bpptgrow2_g) * (pptgrow_seq^2) +
                         mean(grow_coef$bpptdorm2_g) * (grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s])^2 +
                         mean(grow_coef$btempgrow2_g) * (grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s])^2 +
                         mean(grow_coef$btempdorm2_g) * (grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s])^2 +
                         mean(grow_coef$bpptgrow2sex_g) * (pptgrow_seq^2) * (s-1) +
                         mean(grow_coef$bpptdorm2sex_g) * (grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s])^2 * (s-1) +
                         mean(grow_coef$btempgrow2sex_g) * (grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s])^2 * (s-1) +
                         mean(grow_coef$btempdorm2sex_g) * (grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s])^2 * (s-1)
    
    ),lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_grow_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(tempgrow[size_bin==i]*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),mean_grow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n",ylim=c(0,50));box()    
     mylabel3 <- paste(graph[[6]])
     mtext(mylabel3,side = 3, adj = 0,cex=1.5)    
    for(s in 1:2){
      points(tempgrow[Sex==s & size_bin==i]*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),mean_grow[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(tempgrow_seq*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow) ,
            exp(mean(grow_coef$b0_g) +
                        mean(grow_coef$bsize_g) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i] +
                        mean(grow_coef$bsizesex_g) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                        mean(grow_coef$bsex_g) * (s-1) +
                        mean(grow_coef$bpptgrow_g) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s]  +
                        mean(grow_coef$bpptdorm_g) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s]  +
                        mean(grow_coef$btempgrow_g) * tempgrow_seq +
                        mean(grow_coef$btempdorm_g) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] +
                        mean(grow_coef$bpptgrowsex_g) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * (s-1) +
                        mean(grow_coef$bpptdormsex_g) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s]  * (s-1) +
                        mean(grow_coef$btempdormpptdorm_g) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] +
                        mean(grow_coef$btempgrowpptgrow_g) * tempgrow_seq * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] +
                        mean(grow_coef$btempgrowsex_g) * tempgrow_seq * (s-1) +
                        mean(grow_coef$btempdormsex_g) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * (s-1) +
                        mean(grow_coef$btempdormpptdormsex_g) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] * (s-1) +
                        mean(grow_coef$btempgrowpptgrowsex_g) * tempgrow_seq * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s]  * (s-1) +
                        mean(grow_coef$bpptgrow2_g) * (grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s])^2 +
                        mean(grow_coef$bpptdorm2_g) * (grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s])^2 +
                        mean(grow_coef$btempgrow2_g) * (tempgrow_seq^2) +
                        mean(grow_coef$btempdorm2_g) * (grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s])^2 +
                        mean(grow_coef$bpptgrow2sex_g) * (grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s])^2 * (s-1) +
                        mean(grow_coef$bpptdorm2sex_g) * (grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s])^2 * (s-1) +
                        mean(grow_coef$btempgrow2sex_g) * (tempgrow_seq^2) * (s-1) +
                        mean(grow_coef$btempdorm2sex_g) * (grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s])^2 * (s-1)

    
    ),lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})


with(poar_grow_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(pptdorm[size_bin==i]*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),mean_grow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n",ylim=c(0,50));box()    
    # if(i==1){mtext("#tillers",side=2,line=3,cex=1.3)}
    mylabel3 <- paste(graph[[7]])
    mtext(mylabel3,side = 3, adj = 0,cex=1.5)
    for(s in 1:2){
      points(pptdorm[Sex==s & size_bin==i]*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),mean_grow[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2,ylim=c(0,50))
      lines(pptdorm_seq*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm) ,
            exp(mean(grow_coef$b0_g) + 
                         mean(grow_coef$bsize_g) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i]+
                         mean(grow_coef$bsizesex_g) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                         mean(grow_coef$bsex_g) * (s-1) +
                         mean(grow_coef$bpptgrow_g) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] +
                         mean(grow_coef$bpptdorm_g) *  pptdorm_seq +
                         mean(grow_coef$btempgrow_g) * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] +
                         mean(grow_coef$btempdorm_g) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] +
                         mean(grow_coef$bpptgrowsex_g) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * (s-1) +
                         mean(grow_coef$bpptdormsex_g) * pptdorm_seq * (s-1) +
                         mean(grow_coef$btempgrowsex_g) * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] * (s-1) +
                        mean(grow_coef$btempdormsex_g) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * (s-1) +
                        mean (grow_coef$btempdormpptdorm_g) * pptdorm_seq * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] +
                        mean (grow_coef$btempgrowpptgrow_g) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] +
                         mean(grow_coef$btempdormpptdormsex_g) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * pptdorm_seq * (s-1) +
                         mean(grow_coef$btempgrowpptgrowsex_g) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] * (s-1) +
                         mean(grow_coef$bpptgrow2_g) * (grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s])^2 +
                         mean(grow_coef$bpptdorm2_g) * (pptdorm_seq^2) +
                         mean(grow_coef$btempgrow2_g) * (grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s])^2 +
                        mean (grow_coef$btempdorm2_g) * (grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s])^2 +
                         mean(grow_coef$bpptgrow2sex_g) * (grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s])^2 * (s-1) +
                         mean(grow_coef$bpptdorm2sex_g) * (pptdorm_seq^2) * (s-1) +
                         mean(grow_coef$btempgrow2sex_g) * (grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s])^2 * (s-1) +
                         mean(grow_coef$btempdorm2sex_g) * (grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s])^2 * (s-1)
    ),lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_grow_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(tempdorm[size_bin==i]*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mean_grow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n",ylim=c(0,50));box()    
    # if(i==1){mtext("#tillers",side=2,line=3,cex=1.3)}
    mylabel3 <- paste(graph[[8]])
    mtext(mylabel3,side = 3, adj = 0,cex=1.5)   
    for(s in 1:2){
      points(tempdorm[Sex==s & size_bin==i]*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm) ,mean_grow[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2,ylim=c(0,50))
      lines(tempdorm_seq *sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm) ,exp(mean(grow_coef$b0_g) + 
                         mean(grow_coef$bsize_g) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i] +
                         mean(grow_coef$bsizesex_g) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                         mean(grow_coef$bsex_g) * (s-1) +
                         mean(grow_coef$bpptgrow_g) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] +
                        mean(grow_coef$bpptdorm_g) *  grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] +
                         mean(grow_coef$btempgrow_g) * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] +
                         mean(grow_coef$btempdorm_g) * tempdorm_seq +
                         mean(grow_coef$bpptgrowsex_g) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * (s-1) +
                         mean(grow_coef$bpptdormsex_g) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] * (s-1) +
                         mean(grow_coef$btempgrowsex_g) * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] * (s-1) +
                        mean(grow_coef$btempdormsex_g) * tempdorm_seq * (s-1) +
                         mean(grow_coef$btempdormpptdorm_g) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] +
                         mean(grow_coef$btempgrowpptgrow_g) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] +
                         mean(grow_coef$btempdormpptdormsex_g) * tempdorm_seq * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] * (s-1) +
                         mean(grow_coef$btempgrowpptgrowsex_g) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] * (s-1) +
                         mean(grow_coef$bpptgrow2_g) * (grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s])^2 +
                         mean(grow_coef$bpptdorm2_g) * (grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s])^2 +
                         mean(grow_coef$btempgrow2_g) * (grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s])^2 +
                         mean(grow_coef$btempdorm2_g) * (tempdorm_seq^2) +
                         mean(grow_coef$bpptgrow2sex_g) * (grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s])^2 * (s-1) +
                         mean(grow_coef$bpptdorm2sex_g) * (grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s])^2 * (s-1) +
                         mean(grow_coef$btempgrow2sex_g) * (grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s])^2 * (s-1) +
                         mean(grow_coef$btempdorm2sex_g) * (tempdorm_seq^2) * (s-1)
    )
            ,lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_flow_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(pptgrow[size_bin==i]*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),mean_flow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n");box()
    if(i==1){mtext("Pr(flowering)",side=2,line=3,cex=1.5)}
    mylabel3 <- paste(graph[[9]])
    mtext(mylabel3,side = 3, adj = 0,cex=1.5)
    for(s in 1:2){
      points(pptgrow[Sex==s & size_bin==i]*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),mean_flow[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(pptgrow_seq *sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow) ,invlogit(mean(flow_coef$b0_f) + 
                         mean(flow_coef$bsize_f) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i]+
                         mean(flow_coef$bsizesex_f) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                         mean(flow_coef$bsex_f) * (s-1) +
                         mean(flow_coef$bpptgrow_f) * pptgrow_seq +
                         mean(flow_coef$bpptdorm_f) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s]  +
                         mean(flow_coef$btempgrow_f) * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] +
                         mean(flow_coef$btempdorm_f) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] +
                         mean(flow_coef$bpptgrowsex_f) * pptgrow_seq * (s-1) +
                         mean(flow_coef$bpptdormsex_f) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s]  * (s-1) +
                         mean(flow_coef$btempgrowsex_f) * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] * (s-1) +
                        mean(flow_coef$btempdormsex_f) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * (s-1) +
                         mean(flow_coef$btempdormpptdorm_f) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] +
                         mean(flow_coef$btempgrowpptgrow_f) * pptgrow_seq * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] +
                         mean(flow_coef$btempdormpptdormsex_f) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] * (s-1) +
                         mean(flow_coef$btempgrowpptgrowsex_f) * pptgrow_seq * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s]  * (s-1) +
                         mean(flow_coef$bpptgrow2_f) * (pptgrow_seq^2) +
                         mean(flow_coef$bpptdorm2_f) * (flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s])^2 +
                         mean(flow_coef$btempgrow2_f) * (flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s])^2 +
                         mean(flow_coef$btempdorm2_f) * (flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s])^2 +
                         mean(flow_coef$bpptgrow2sex_f) * (pptgrow_seq^2) * (s-1) +
                         mean(flow_coef$bpptdorm2sex_f) * (flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s])^2 * (s-1) +
                         mean(flow_coef$btempgrow2sex_f) * (flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s])^2 * (s-1) +
                         mean(flow_coef$btempdorm2sex_f) * (flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s])^2 * (s-1)
    
    )
    
            ,lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_flow_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(tempgrow[size_bin==i]*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow) ,mean_flow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n");box()
    mylabel3 <- paste(graph[[10]])
    mtext(mylabel3,side = 3, adj = 0,cex=1.5)
    for(s in 1:2){
      points(tempgrow[Sex==s & size_bin==i]*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),mean_flow[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(tempgrow_seq *sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),invlogit(mean(flow_coef$b0_f) +
                        mean(flow_coef$bsize_f) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] +
                        mean(flow_coef$bsizesex_f) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                        mean(flow_coef$bsex_f) * (s-1) +
                        mean(flow_coef$bpptgrow_f) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s]  +
                        mean(flow_coef$bpptdorm_f) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s]  +
                        mean(flow_coef$btempgrow_f) * tempgrow_seq +
                        mean(flow_coef$btempdorm_f) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] +
                        mean(flow_coef$bpptgrowsex_f) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * (s-1) +
                        mean(flow_coef$bpptdormsex_f) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s]  * (s-1) +
                        mean(flow_coef$btempgrowsex_f) * tempgrow_seq * (s-1) +
                        mean(flow_coef$btempdormsex_f) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * (s-1) +
                        mean(flow_coef$btempdormpptdorm_f) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s]  +
                        mean(flow_coef$btempgrowpptgrow_f) * tempgrow_seq * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s]   +
                        mean(flow_coef$btempdormpptdormsex_f) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] * (s-1) +
                        mean(flow_coef$btempgrowpptgrowsex_f) * tempgrow_seq * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s]  * (s-1) +
                        mean(flow_coef$bpptgrow2_f) * (flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s])^2 +
                        mean(flow_coef$bpptdorm2_f) * (flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s])^2 +
                        mean(flow_coef$btempgrow2_f) * (tempgrow_seq^2) +
                        mean(flow_coef$btempdorm2_f) * (flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s])^2 +
                        mean(flow_coef$bpptgrow2sex_f) * (flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s])^2 * (s-1) +
                        mean(flow_coef$bpptdorm2sex_f) * (flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s])^2 * (s-1) +
                        mean(flow_coef$btempgrow2sex_f) * (tempgrow_seq^2) * (s-1) +
                        mean(flow_coef$btempdorm2sex_f) * (flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s])^2 * (s-1)
    ) ,lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_flow_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(pptdorm[size_bin==i]*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm) ,mean_flow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n");box()
     mylabel3 <- paste(graph[[11]])
     mtext(mylabel3,side = 3, adj = 0,cex=1.5)  
    for(s in 1:2){
      points(pptdorm[Sex==s & size_bin==i]*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),mean_flow[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(pptdorm_seq*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm), invlogit(mean(flow_coef$b0_f) + 
                         mean(flow_coef$bsize_f) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i]+
                         mean(flow_coef$bsizesex_f) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                         mean(flow_coef$bsex_f) * (s-1) +
                         mean(flow_coef$bpptgrow_f) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] +
                         mean(flow_coef$bpptdorm_f) *  pptdorm_seq +
                         mean(flow_coef$btempgrow_f) * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] +
                         mean(flow_coef$btempdorm_f) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] +
                         mean(flow_coef$bpptgrowsex_f) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * (s-1) +
                         mean(flow_coef$bpptdormsex_f) * pptdorm_seq * (s-1) +
                        mean (flow_coef$btempgrowsex_f) * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] * (s-1) +
                        mean(flow_coef$btempdormsex_f) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * (s-1) +
                        mean (flow_coef$btempdormpptdorm_f) * pptdorm_seq * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] +
                         mean(flow_coef$btempgrowpptgrow_f) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] +
                         mean(flow_coef$btempdormpptdormsex_f) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * pptdorm_seq * (s-1) +
                        mean (flow_coef$btempgrowpptgrowsex_f) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] * (s-1) +
                         mean(flow_coef$bpptgrow2_f) * (flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s])^2 +
                         mean(flow_coef$bpptdorm2_f) * (pptdorm_seq^2) +
                         mean(flow_coef$btempgrow2_f) * (flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s])^2 +
                        mean (flow_coef$btempdorm2_f) * (flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s])^2 +
                         mean(flow_coef$bpptgrow2sex_f) * (flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s])^2 * (s-1) +
                         mean(flow_coef$bpptdorm2sex_f) * (pptdorm_seq^2) * (s-1) +
                        mean (flow_coef$btempgrow2sex_f) * (flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s])^2 * (s-1) +
                        mean (flow_coef$btempdorm2sex_f) * (flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s])^2 * (s-1)
    )
  
    
            ,lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_flow_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(0,4,2,0))
    plot(tempdorm[size_bin==i]*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm) ,mean_flow[size_bin==i],type="n",
         xlab=" ",ylab=" ",xaxt="n");box()
    # if(i==1){mtext("Pr(flowering)",side=2,line=3,cex=1.3)}
     mylabel3 <- paste(graph[[12]])
    mtext(mylabel3,side = 3, adj = 0,cex=1.5)    
    for(s in 1:2){
      points(tempdorm[Sex==s & size_bin==i]*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mean_flow[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=5*(bin_n/max(bin_n)),lwd=2)
      lines(tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm), invlogit(mean(flow_coef$b0_f) + 
                         mean(flow_coef$bsize_f) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] +
                         mean(flow_coef$bsizesex_f) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                         mean(flow_coef$bsex_f) * (s-1) +
                         mean(flow_coef$bpptgrow_f) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] +
                         mean(flow_coef$bpptdorm_f) *  flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] +
                         mean(flow_coef$btempgrow_f) * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] +
                         mean(flow_coef$btempdorm_f) * tempdorm_seq +
                         mean(flow_coef$bpptgrowsex_f) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * (s-1) +
                         mean(flow_coef$bpptdormsex_f) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] * (s-1) +
                        mean (flow_coef$btempgrowsex_f) * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] * (s-1) +
                        mean(flow_coef$btempdormsex_f) * tempdorm_seq * (s-1) +
                        mean (flow_coef$btempdormpptdorm_f) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] +
                         mean(flow_coef$btempgrowpptgrow_f) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] +
                        mean (flow_coef$btempdormpptdormsex_f) * tempdorm_seq * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] * (s-1) +
                        mean (flow_coef$btempgrowpptgrowsex_f) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] * (s-1) +
                        mean (flow_coef$bpptgrow2_f) * (flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s])^2 +
                        mean (flow_coef$bpptdorm2_f) * (flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s])^2 +
                        mean (flow_coef$btempgrow2_f) * (flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s])^2 +
                        mean (flow_coef$btempdorm2_f) * (tempdorm_seq^2) +
                        mean (flow_coef$bpptgrow2sex_f) * (flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s])^2 * (s-1) +
                        mean (flow_coef$bpptdorm2sex_f) * (flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s])^2 * (s-1) +
                        mean (flow_coef$btempgrow2sex_f) * (flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s])^2 * (s-1) +
                        mean (flow_coef$btempdorm2sex_f) * (tempdorm_seq^2) * (s-1)
    )
  
    
            ,lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_panic_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(1,4,2,0))
    plot(pptgrow[size_bin==i]*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow) ,mean_panic[size_bin==i],type="n",
         xlab="",ylab=" ",ylim=c(0,15));box()    
    if(i==1){mtext("#panicles",side=2,line=3,cex=1.5)}
    mylabel3 <- paste(graph[[13]])
    mtext(mylabel3,side = 3, adj = 0,cex=1.75) 
    mtext("Growing season precip",side=1,line=3,cex=1.3)
    for(s in 1:2){
      points(pptgrow[Sex==s & size_bin==i]*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow),mean_panic[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=7*(bin_n/max(bin_n)),lwd=2,ylim=c(0,12))
      lines(pptgrow_seq*sd(poar_2015_2016$pptgrow) + mean(poar_2015_2016$pptgrow), exp(mean(panic_coef$b0_p) + 
                         mean(panic_coef$bsize_p) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i]+
                         mean(panic_coef$bsizesex_p) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                         mean(panic_coef$bsex_p) * (s-1) +
                         mean(panic_coef$bpptgrow_p) * pptgrow_seq +
                         mean(panic_coef$bpptdorm_p) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s]  +
                         mean(panic_coef$btempgrow_p) * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] +
                         mean(panic_coef$btempdorm_p) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] +
                         mean(panic_coef$bpptgrowsex_p) * pptgrow_seq * (s-1) +
                         mean(panic_coef$bpptdormsex_p) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s]  * (s-1) +
                         mean(panic_coef$btempgrowsex_p) * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] * (s-1) +
                        mean(panic_coef$btempdormsex_p) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * (s-1) +
                         mean(panic_coef$btempdormpptdorm_p) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] +
                         mean(panic_coef$btempgrowpptgrow_p) * pptgrow_seq * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] +
                         mean(panic_coef$btempdormpptdormsex_p) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] * (s-1) +
                         mean(panic_coef$btempgrowpptgrowsex_p) * pptgrow_seq * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s]  * (s-1) +
                         mean(panic_coef$bpptgrow2_p) * (pptgrow_seq^2) +
                         mean(panic_coef$bpptdorm2_p) * (panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s])^2 +
                         mean(panic_coef$btempgrow2_p) * (panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s])^2 +
                         mean(panic_coef$btempdorm2_p) * (panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s])^2 +
                         mean(panic_coef$bpptgrow2sex_p) * (pptgrow_seq^2) * (s-1) +
                         mean(panic_coef$bpptdorm2sex_p) * (panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s])^2 * (s-1) +
                         mean(panic_coef$btempgrow2sex_p) * (panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s])^2 * (s-1) +
                         mean(panic_coef$btempdorm2sex_p) * (panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s])^2 * (s-1)
    
    )
    
            ,lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
    legend(900, 12, 
       legend=c( "Female","Male"),
       lty=c(1,1),
       # pt.cex=c(0.75,1.25,1),
       col =c("#D95F02", "#1B9E77"),
       cex = 1.2,
       bty = "n",
       lw=3,
       horiz = F ) 
  }
})

with(poar_panic_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(1,4,2,0))
    plot(tempgrow[size_bin==i]* sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),mean_panic[size_bin==i],type="n",
         xlab=" ",ylab=" ",ylim=c(0,12));box()    
    # if(i==1){mtext("#panicles",side=2,line=3,cex=1.3)}
    mylabel3 <- paste(graph[[14]])
    mtext(mylabel3,side = 3, adj = 0,cex=1.5)    
    mtext("Growing season temp",side=1,line=3,cex=1.3)
    for(s in 1:2){
      points(tempgrow[Sex==s & size_bin==i]*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow),mean_panic[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=7*(bin_n/max(bin_n)),lwd=2)
      lines(tempgrow_seq*sd(poar_2015_2016$tempgrow) + mean(poar_2015_2016$tempgrow), exp(mean(panic_coef$b0_p) +
                        mean(panic_coef$bsize_p) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i] +
                        mean(panic_coef$bsizesex_p) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                        mean(panic_coef$bsex_p) * (s-1) +
                        mean(panic_coef$bpptgrow_p) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s]  +
                        mean(panic_coef$bpptdorm_p) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s]  +
                        mean(panic_coef$btempgrow_p) * tempgrow_seq +
                        mean(panic_coef$btempdorm_p) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] +
                        mean(panic_coef$bpptgrowsex_p) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] * (s-1) +
                        mean(panic_coef$bpptdormsex_p) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s]  * (s-1) +
                        mean(panic_coef$btempgrowsex_p) * tempgrow_seq * (s-1) +
                        mean(panic_coef$btempdormsex_p) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * (s-1) +
                        mean(panic_coef$btempdormpptdorm_p) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s]  +
                        mean(panic_coef$btempgrowpptgrow_p) * tempgrow_seq * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s]   +
                        mean(panic_coef$btempdormpptdormsex_p) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] * (s-1) +
                        mean(panic_coef$btempgrowpptgrowsex_p) * tempgrow_seq * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s]  * (s-1) +
                        mean(panic_coef$bpptgrow2_p) * (panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s])^2 +
                        mean(panic_coef$bpptdorm2_p) * (panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s])^2 +
                        mean(panic_coef$btempgrow2_p) * (tempgrow_seq^2) +
                        mean(panic_coef$btempdorm2_p) * (panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s])^2 +
                        mean(panic_coef$bpptgrow2sex_p) * (panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s])^2 * (s-1) +
                        mean(panic_coef$bpptdorm2sex_p) * (panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s])^2 * (s-1) +
                        mean(panic_coef$btempgrow2sex_p) * (tempgrow_seq^2) * (s-1) +
                        mean(panic_coef$btempdorm2sex_p) * (panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s])^2 * (s-1)
    
    ) ,lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_panic_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(1,4,2,0))
    plot(pptdorm[size_bin==i]*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),mean_panic[size_bin==i],type="n",
         xlab=" ",ylab=" ",ylim=c(0,12));box()    
    # if(i==1){mtext("#panicles",side=2,line=3,cex=1.3)}
     mylabel3 <- paste(graph[[15]])
     mtext(mylabel3,side = 3, adj = 0,cex=1.5)    
     mtext("Dormant season precip",side=1,line=3,cex=1.3)
    for(s in 1:2){
      points(pptdorm[Sex==s & size_bin==i]*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm),mean_panic[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=7*(bin_n/max(bin_n)),lwd=2)
      lines(pptdorm_seq*sd(poar_2015_2016$pptdorm) + mean(poar_2015_2016$pptdorm), exp(mean(panic_coef$b0_p) + 
                         mean(panic_coef$bsize_p) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i]+
                         mean(panic_coef$bsizesex_p) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                         mean(panic_coef$bsex_p) * (s-1) +
                         mean(panic_coef$bpptdorm_p) * pptdorm_seq +
                         mean(panic_coef$bpptgrow_p) * panic_mean_pptgrow$pptgrow[panic_mean_pptdorm$Sex==s]  +
                         mean(panic_coef$btempgrow_p) * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] +
                         mean(panic_coef$btempdorm_p) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] +
                         mean(panic_coef$bpptdormsex_p) * pptdorm_seq * (s-1) +
                         mean(panic_coef$bpptgrowsex_p) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s]  * (s-1) +
                         mean(panic_coef$btempgrowsex_p) * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] * (s-1) +
                        mean(panic_coef$btempdormsex_p) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * (s-1) +
                         mean(panic_coef$btempdormpptdorm_p) * pptdorm_seq * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] +
                         mean(panic_coef$btempgrowpptgrow_p) * panic_mean_pptgrow$pptgrow[panic_mean_pptdorm$Sex==s] * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] +
                         mean(panic_coef$btempdormpptdormsex_p) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * pptdorm_seq * (s-1) +
                         mean(panic_coef$btempgrowpptgrowsex_p) * panic_mean_pptgrow$pptgrow[panic_mean_tempgrow$Sex==s] * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s]  * (s-1) +
                         mean(panic_coef$bpptdorm2_p) * (pptdorm_seq^2) +
                         mean(panic_coef$bpptgrow2_p) * (panic_mean_pptgrow$pptgrow[panic_mean_pptdorm$Sex==s])^2 +
                         mean(panic_coef$btempgrow2_p) * (panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s])^2 +
                         mean(panic_coef$btempdorm2_p) * (panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s])^2 +
                         mean(panic_coef$bpptdorm2sex_p) * (pptdorm_seq^2) * (s-1) +
                         mean(panic_coef$bpptgrow2sex_p) * (panic_mean_pptgrow$pptgrow[panic_mean_tempgrow$Sex==s])^2 * (s-1) +
                         mean(panic_coef$btempgrow2sex_p) * (panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s])^2 * (s-1) +
                         mean(panic_coef$btempdorm2sex_p) * (panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s])^2 * (s-1)
    
    )
    
            ,lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

with(poar_panic_binned,{
  for(i in 1:size_bin_num){
    par(mar=c(1,4,2,0))
    plot(tempdorm[size_bin==i]*sd(poar_2015_2016$tempdorm)+ mean(poar_2015_2016$tempdorm) ,mean_panic[size_bin==i],type="n",
         xlab=" ",ylab=" ",ylim=c(0,12));box()    
    # if(i==1){mtext("#panicles",side=2,line=3,cex=1.3)}
    mylabel3 <- paste(graph[[16]])
    mtext(mylabel3,side = 3, adj = 0,cex=1.5)   
    mtext("Dormant season temp",side=1,line=3,cex=1.3)
    for(s in 1:2){
      points(tempdorm[Sex==s & size_bin==i]*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mean_panic[Sex==s & size_bin==i],
             bg=sex_cols[s],pch=21,cex=7*(bin_n/max(bin_n)),lwd=2)
      lines(tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm), exp(mean(panic_coef$b0_p) + 
                         mean(panic_coef$bsize_p) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i] +
                         mean(panic_coef$bsizesex_p) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                         mean(panic_coef$bsex_p) * (s-1) +
                         mean(panic_coef$bpptgrow_p) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] +
                         mean(panic_coef$bpptdorm_p) *  panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] +
                         mean(panic_coef$btempgrow_p) * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] +
                         mean(panic_coef$btempdorm_p) * tempdorm_seq +
                         mean(panic_coef$bpptgrowsex_p) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] * (s-1) +
                         mean(panic_coef$bpptdormsex_p) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] * (s-1) +
                         mean(panic_coef$btempgrowsex_p) * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] * (s-1) +
                        mean(panic_coef$btempdormsex_p) * tempdorm_seq * (s-1) +
                        mean (panic_coef$btempdormpptdorm_p) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] * tempdorm_seq +
                        mean (panic_coef$btempgrowpptgrow_p) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] +
                        mean (panic_coef$btempdormpptdormsex_p) * tempdorm_seq * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] * (s-1) +
                        mean (panic_coef$btempgrowpptgrowsex_p) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] * (s-1) +
                        mean (panic_coef$bpptgrow2_p) * (panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s])^2 +
                        mean (panic_coef$bpptdorm2_p) * (panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s])^2 +
                         mean(panic_coef$btempgrow2_p) * (panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s])^2 +
                         mean(panic_coef$btempdorm2_p) * (tempdorm_seq^2) +
                         mean(panic_coef$bpptgrow2sex_p) * (panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s])^2 * (s-1) +
                         mean(panic_coef$bpptdorm2sex_p) * (panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s])^2 * (s-1) +
                         mean(panic_coef$btempgrow2sex_p) * (panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s])^2 * (s-1) +
                         mean(panic_coef$btempdorm2sex_p) * (tempdorm_seq^2) * (s-1)
    )
     ,lty=sex_lty[s],lwd=3,col= sex_cols[s])
    }
  }
})

dev.off()

# Plot vital rate figure (3D) ----
## Growing season---- 
pptgrow_seq <- seq(min(poar_surv_binned$pptgrow),max(poar_surv_binned$pptgrow),length.out=30)
tempgrow_seq <- seq(min(poar_surv_binned$tempgrow),max(poar_surv_binned$tempgrow),length.out=30)
climgrow<-expand.grid(ppgrow=pptgrow_seq,tempgrow=tempgrow_seq)
n_post_draws <- 2000
post_draws <- sample.int(length(surv_coef$b0_s), n_post_draws)
surv_sex_diff_post_grow <- grow_sex_diff_post_grow <- flow_sex_diff_post_grow <- panic_sex_diff_post_grow <- array(NA,dim=c(size_bin_num,nrow(climgrow),n_post_draws))
for(p in 1:n_post_draws){
  for(i in 1:size_bin_num){
    s=1;      
    fem_s <- invlogit((surv_coef$b0_s[post_draws[p]]) + 
                        (surv_coef$bsize_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i]+
                        (surv_coef$bsizesex_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                        (surv_coef$bsex_s[post_draws[p]]) * (s-1) +
                        (surv_coef$bpptgrow_s[post_draws[p]]) * climgrow$ppgrow +
                        (surv_coef$bpptdorm_s[post_draws[p]]) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s]  +
                        (surv_coef$btempgrow_s[post_draws[p]]) * climgrow$tempgrow +
                        (surv_coef$btempdorm_s[post_draws[p]]) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] +
                        (surv_coef$bpptgrowsex_s[post_draws[p]]) * climgrow$ppgrow * (s-1) +
                        (surv_coef$bpptdormsex_s[post_draws[p]]) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s]  * (s-1) +
                        (surv_coef$btempgrowsex_s[post_draws[p]]) * climgrow$tempgrow * (s-1) +
                        (surv_coef$btempdormsex_s[post_draws[p]]) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * (s-1) +
                        (surv_coef$btempdormpptdorm_s[post_draws[p]]) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] +
                        (surv_coef$btempgrowpptgrow_s[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow +
                        (surv_coef$btempdormpptdormsex_s[post_draws[p]]) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] * (s-1) +
                        (surv_coef$btempgrowpptgrowsex_s[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow * (s-1) +
                        (surv_coef$bpptgrow2_s[post_draws[p]]) * (climgrow$ppgrow)^2 +
                        (surv_coef$bpptdorm2_s[post_draws[p]]) * (surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s])^2 +
                        (surv_coef$btempgrow2_s[post_draws[p]]) * (climgrow$tempgrow)^2 +
                        (surv_coef$btempdorm2_s[post_draws[p]]) * (surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s])^2 +
                        (surv_coef$bpptgrow2sex_s[post_draws[p]]) * (climgrow$ppgrow)^2 * (s-1) +
                        (surv_coef$bpptdorm2sex_s[post_draws[p]]) * (surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s])^2 * (s-1) +
                        (surv_coef$btempgrow2sex_s[post_draws[p]]) * (climgrow$tempgrow)^2 * (s-1) +
                        (surv_coef$btempdorm2sex_s[post_draws[p]]) * (surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s])^2 * (s-1)
    )
    
    fem_g <- exp((grow_coef$b0_g[post_draws[p]]) + 
                   (grow_coef$bsize_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i]+
                   (grow_coef$bsizesex_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                   (grow_coef$bsex_g[post_draws[p]]) * (s-1) +
                   (grow_coef$bpptgrow_g[post_draws[p]]) * climgrow$ppgrow +
                   (grow_coef$bpptdorm_g[post_draws[p]]) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s]  +
                   (grow_coef$btempgrow_g[post_draws[p]]) * climgrow$tempgrow +
                   (grow_coef$btempdorm_g[post_draws[p]]) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] +
                   (grow_coef$bpptgrowsex_g[post_draws[p]]) * climgrow$ppgrow * (s-1) +
                   (grow_coef$bpptdormsex_g[post_draws[p]]) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s]  * (s-1) +
                   (grow_coef$btempgrowsex_g[post_draws[p]]) * climgrow$tempgrow * (s-1) +
                   (grow_coef$btempdormsex_g[post_draws[p]]) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * (s-1) +
                   (grow_coef$btempdormpptdorm_g[post_draws[p]]) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] +
                   (grow_coef$btempgrowpptgrow_g[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow +
                   (grow_coef$btempdormpptdormsex_g[post_draws[p]]) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] * (s-1) +
                   (grow_coef$btempgrowpptgrowsex_g[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow  * (s-1) +
                   (grow_coef$bpptgrow2_g[post_draws[p]]) * (climgrow$ppgrow)^2 +
                   (grow_coef$bpptdorm2_g[post_draws[p]]) * (grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s])^2 +
                   (grow_coef$btempgrow2_g[post_draws[p]]) * (climgrow$tempgrow)^2 +
                   (grow_coef$btempdorm2_g[post_draws[p]]) * (grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s])^2 +
                   (grow_coef$bpptgrow2sex_g[post_draws[p]]) * (climgrow$ppgrow)^2 * (s-1) +
                   (grow_coef$bpptdorm2sex_g[post_draws[p]]) * (grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s])^2 * (s-1) +
                   (grow_coef$btempgrow2sex_g[post_draws[p]]) * (climgrow$tempgrow)^2 * (s-1) +
                   (grow_coef$btempdorm2sex_g[post_draws[p]]) * (grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s])^2 * (s-1)
    )
    
    fem_f <- invlogit((flow_coef$b0_f[post_draws[p]]) + 
                        (flow_coef$bsize_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] +
                        (flow_coef$bsizesex_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                        (flow_coef$bsex_f[post_draws[p]]) * (s-1) +
                        (flow_coef$bpptgrow_f[post_draws[p]]) * climgrow$ppgrow +
                        (flow_coef$bpptdorm_f[post_draws[p]]) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s]  +
                        (flow_coef$btempgrow_f[post_draws[p]]) * climgrow$tempgrow +
                        (flow_coef$btempdorm_f[post_draws[p]]) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] +
                        (flow_coef$bpptgrowsex_f[post_draws[p]]) * climgrow$ppgrow * (s-1) +
                        (flow_coef$bpptdormsex_f[post_draws[p]]) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s]  * (s-1) +
                        (flow_coef$btempgrowsex_f[post_draws[p]]) * climgrow$tempgrow * (s-1) +
                        (flow_coef$btempdormsex_f[post_draws[p]]) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * (s-1) +
                        (flow_coef$btempdormpptdorm_f[post_draws[p]]) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] +
                        (flow_coef$btempgrowpptgrow_f[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow +
                        (flow_coef$btempdormpptdormsex_f[post_draws[p]]) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] * (s-1) +
                        (flow_coef$btempgrowpptgrowsex_f)[post_draws[p]] * climgrow$ppgrow * climgrow$tempgrow  * (s-1) +
                        (flow_coef$bpptgrow2_f[post_draws[p]]) * (climgrow$ppgrow)^2 +
                        (flow_coef$bpptdorm2_f[post_draws[p]]) * (flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s])^2 +
                        (flow_coef$btempgrow2_f[post_draws[p]]) * (climgrow$tempgrow)^2 +
                        (flow_coef$btempdorm2_f[post_draws[p]]) * (flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s])^2 +
                        (flow_coef$bpptgrow2sex_f[post_draws[p]]) * (climgrow$ppgrow)^2 * (s-1) +
                        (flow_coef$bpptdorm2sex_f[post_draws[p]]) * (flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s])^2 * (s-1) +
                        (flow_coef$btempgrow2sex_f[post_draws[p]]) * (climgrow$tempgrow)^2 * (s-1) +
                        (flow_coef$btempdorm2sex_f[post_draws[p]]) * (flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s])^2 * (s-1)
                      
    )
    
    fem_p <- exp((panic_coef$b0_p[post_draws[p]]) + 
                   (panic_coef$bsize_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i]+
                   (panic_coef$bsizesex_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                   (panic_coef$bsex_p[post_draws[p]]) * (s-1) +
                   (panic_coef$bpptgrow_p[post_draws[p]]) * climgrow$ppgrow +
                   (panic_coef$bpptdorm_p[post_draws[p]]) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s]  +
                   (panic_coef$btempgrow_p[post_draws[p]]) * climgrow$tempgrow  +
                   (panic_coef$btempdorm_p[post_draws[p]]) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] +
                   (panic_coef$bpptgrowsex_p[post_draws[p]]) * climgrow$ppgrow * (s-1) +
                   (panic_coef$bpptdormsex_p[post_draws[p]]) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s]  * (s-1) +
                   (panic_coef$btempgrowsex_p[post_draws[p]]) * climgrow$tempgrow  * (s-1) +
                   (panic_coef$btempdormsex_p[post_draws[p]]) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * (s-1) +
                   (panic_coef$btempdormpptdorm_p[post_draws[p]]) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] +
                   (panic_coef$btempgrowpptgrow_p[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow  +
                   (panic_coef$btempdormpptdormsex_p[post_draws[p]]) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] * (s-1) +
                   (panic_coef$btempgrowpptgrowsex_p[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow  * (s-1) +
                   (panic_coef$bpptgrow2_p[post_draws[p]]) * (climgrow$ppgrow)^2 +
                   (panic_coef$bpptdorm2_p[post_draws[p]]) * (panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s])^2 +
                   (panic_coef$btempgrow2_p[post_draws[p]]) * (climgrow$tempgrow)^2 +
                   (panic_coef$btempdorm2_p[post_draws[p]]) * (panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s])^2 +
                   (panic_coef$bpptgrow2sex_p[post_draws[p]]) * (climgrow$ppgrow)^2 * (s-1) +
                   (panic_coef$bpptdorm2sex_p[post_draws[p]]) * (panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s])^2 * (s-1) +
                   (panic_coef$btempgrow2sex_p[post_draws[p]]) * (climgrow$tempgrow)^2 * (s-1) +
                   (panic_coef$btempdorm2sex_p[post_draws[p]]) * (panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s])^2 * (s-1)
    )
    
    s=2;       
    male_s <- invlogit((surv_coef$b0_s[post_draws[p]]) + 
                         (surv_coef$bsize_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i]+
                         (surv_coef$bsizesex_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                         (surv_coef$bsex_s[post_draws[p]]) * (s-1) +
                         (surv_coef$bpptgrow_s[post_draws[p]]) * climgrow$ppgrow +
                         (surv_coef$bpptdorm_s[post_draws[p]]) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s]  +
                         (surv_coef$btempgrow_s[post_draws[p]]) * climgrow$tempgrow  +
                         (surv_coef$btempdorm_s[post_draws[p]]) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] +
                         (surv_coef$bpptgrowsex_s[post_draws[p]]) * climgrow$ppgrow * (s-1) +
                         (surv_coef$bpptdormsex_s[post_draws[p]]) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s]  * (s-1) +
                         (surv_coef$btempgrowsex_s[post_draws[p]]) * climgrow$tempgrow * (s-1) +
                         (surv_coef$btempdormsex_s[post_draws[p]]) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * (s-1) +
                         (surv_coef$btempdormpptdorm_s[post_draws[p]]) * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] +
                         (surv_coef$btempgrowpptgrow_s[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow  +
                         (surv_coef$btempdormpptdormsex_s[post_draws[p]]) * surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s] * surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s] * (s-1) +
                         (surv_coef$btempgrowpptgrowsex_s[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow   * (s-1) +
                         (surv_coef$bpptgrow2_s[post_draws[p]]) * (climgrow$ppgrow)^2 +
                         (surv_coef$bpptdorm2_s[post_draws[p]]) * (surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s])^2 +
                         (surv_coef$btempgrow2_s[post_draws[p]]) * (climgrow$tempgrow)^2 +
                         (surv_coef$btempdorm2_s[post_draws[p]]) * (surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s])^2 +
                         (surv_coef$bpptgrow2sex_s[post_draws[p]]) * (climgrow$ppgrow)^2 * (s-1) +
                         (surv_coef$bpptdorm2sex_s[post_draws[p]]) * (surv_mean_pptdorm$pptdorm[surv_mean_pptdorm$Sex==s])^2 * (s-1) +
                         (surv_coef$btempgrow2sex_s[post_draws[p]]) * (climgrow$tempgrow )^2 * (s-1) +
                         (surv_coef$btempdorm2sex_s[post_draws[p]]) * (surv_mean_tempdorm$tempdorm[surv_mean_tempdorm$Sex==s])^2 * (s-1)
    )
    
    male_g <- exp((grow_coef$b0_g[post_draws[p]]) + 
                    (grow_coef$bsize_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i]+
                    (grow_coef$bsizesex_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                    (grow_coef$bsex_g[post_draws[p]]) * (s-1) +
                    (grow_coef$bpptgrow_g[post_draws[p]]) * climgrow$ppgrow +
                    (grow_coef$bpptdorm_g[post_draws[p]]) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s]  +
                    (grow_coef$btempgrow_g[post_draws[p]]) * climgrow$tempgrow +
                    (grow_coef$btempdorm_g[post_draws[p]]) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] +
                    (grow_coef$bpptgrowsex_g[post_draws[p]]) * climgrow$ppgrow * (s-1) +
                    (grow_coef$bpptdormsex_g[post_draws[p]]) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s]  * (s-1) +
                    (grow_coef$btempgrowsex_g[post_draws[p]]) * climgrow$tempgrow * (s-1) +
                    (grow_coef$btempdormsex_g[post_draws[p]]) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * (s-1) +
                    (grow_coef$btempdormpptdorm_g[post_draws[p]]) * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] +
                    (grow_coef$btempgrowpptgrow_g[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow +
                    (grow_coef$btempdormpptdormsex_g[post_draws[p]]) * grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s] * grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s] * (s-1) +
                    (grow_coef$btempgrowpptgrowsex_g)[post_draws[p]] * climgrow$ppgrow * climgrow$tempgrow  * (s-1) +
                    (grow_coef$bpptgrow2_g[post_draws[p]]) * (climgrow$ppgrow)^2 +
                    (grow_coef$bpptdorm2_g[post_draws[p]]) * (grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s])^2 +
                    (grow_coef$btempgrow2_g[post_draws[p]]) * (climgrow$tempgrow)^2 +
                    (grow_coef$btempdorm2_g[post_draws[p]]) * (grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s])^2 +
                    (grow_coef$bpptgrow2sex_g[post_draws[p]]) * (climgrow$ppgrow)^2 * (s-1) +
                    (grow_coef$bpptdorm2sex_g[post_draws[p]]) * (grow_mean_pptdorm$pptdorm[grow_mean_pptdorm$Sex==s])^2 * (s-1) +
                    (grow_coef$btempgrow2sex_g[post_draws[p]]) * (climgrow$tempgrow)^2 * (s-1) +
                    (grow_coef$btempdorm2sex_g[post_draws[p]]) * (grow_mean_tempdorm$tempdorm[grow_mean_tempdorm$Sex==s])^2 * (s-1)
                  
    )
    
    male_f <- invlogit((flow_coef$b0_f[post_draws[p]]) + 
                         (flow_coef$bsize_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i]+
                         (flow_coef$bsizesex_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                         (flow_coef$bsex_f[post_draws[p]]) * (s-1) +
                         (flow_coef$bpptgrow_f[post_draws[p]]) * climgrow$ppgrow +
                         (flow_coef$bpptdorm_f[post_draws[p]]) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s]  +
                         (flow_coef$btempgrow_f[post_draws[p]]) * climgrow$tempgrow +
                         (flow_coef$btempdorm_f[post_draws[p]]) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] +
                         (flow_coef$bpptgrowsex_f[post_draws[p]]) * climgrow$ppgrow * (s-1) +
                         (flow_coef$bpptdormsex_f[post_draws[p]]) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s]  * (s-1) +
                         (flow_coef$btempgrowsex_f[post_draws[p]]) * climgrow$tempgrow * (s-1) +
                         (flow_coef$btempdormsex_f[post_draws[p]]) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * (s-1) +
                         (flow_coef$btempdormpptdorm_f[post_draws[p]]) * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] +
                         (flow_coef$btempgrowpptgrow_f[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow +
                         (flow_coef$btempdormpptdormsex_f[post_draws[p]]) * flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s] * flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s] * (s-1) +
                         (flow_coef$btempgrowpptgrowsex_f)[post_draws[p]] * climgrow$ppgrow * climgrow$tempgrow  * (s-1) +
                         (flow_coef$bpptgrow2_f[post_draws[p]]) * (climgrow$ppgrow)^2 +
                         (flow_coef$bpptdorm2_f[post_draws[p]]) * (flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s])^2 +
                         (flow_coef$btempgrow2_f[post_draws[p]]) * (climgrow$tempgrow)^2 +
                         (flow_coef$btempdorm2_f[post_draws[p]]) * (flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s])^2 +
                         (flow_coef$bpptgrow2sex_f[post_draws[p]]) * (climgrow$ppgrow)^2 * (s-1) +
                         (flow_coef$bpptdorm2sex_f[post_draws[p]]) * (flow_mean_pptdorm$pptdorm[flow_mean_pptdorm$Sex==s])^2 * (s-1) +
                         (flow_coef$btempgrow2sex_f[post_draws[p]]) * (climgrow$tempgrow)^2 * (s-1) +
                         (flow_coef$btempdorm2sex_f[post_draws[p]]) * (flow_mean_tempdorm$tempdorm[flow_mean_tempdorm$Sex==s])^2 * (s-1)
    )
    
    male_p <- exp((panic_coef$b0_p[post_draws[p]]) + 
                    (panic_coef$bsize_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i]+
                    (panic_coef$bsizesex_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                    (panic_coef$bsex_p[post_draws[p]]) * (s-1) +
                    (panic_coef$bpptgrow_p[post_draws[p]]) * climgrow$ppgrow +
                    (panic_coef$bpptdorm_p[post_draws[p]]) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s]  +
                    (panic_coef$btempgrow_p[post_draws[p]]) * climgrow$tempgrow +
                    (panic_coef$btempdorm_p[post_draws[p]]) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] +
                    (panic_coef$bpptgrowsex_p[post_draws[p]]) * climgrow$ppgrow * (s-1) +
                    (panic_coef$bpptdormsex_p[post_draws[p]]) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s]  * (s-1) +
                    (panic_coef$btempgrowsex_p[post_draws[p]]) * climgrow$tempgrow * (s-1) +
                    (panic_coef$btempdormsex_p[post_draws[p]]) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * (s-1) +
                    (panic_coef$btempdormpptdorm_p[post_draws[p]]) * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] +
                    (panic_coef$btempgrowpptgrow_p[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow +
                    (panic_coef$btempdormpptdormsex_p[post_draws[p]]) * panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s] * panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s] * (s-1) +
                    (panic_coef$btempgrowpptgrowsex_p[post_draws[p]]) * climgrow$ppgrow * climgrow$tempgrow  * (s-1) +
                    (panic_coef$bpptgrow2_p[post_draws[p]]) * (climgrow$ppgrow)^2 +
                    (panic_coef$bpptdorm2_p[post_draws[p]]) * (panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s])^2 +
                    (panic_coef$btempgrow2_p[post_draws[p]]) * (climgrow$tempgrow)^2 +
                    (panic_coef$btempdorm2_p[post_draws[p]]) * (panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s])^2 +
                    (panic_coef$bpptgrow2sex_p[post_draws[p]]) * (climgrow$ppgrow)^2 * (s-1) +
                    (panic_coef$bpptdorm2sex_p[post_draws[p]]) * (panic_mean_pptdorm$pptdorm[panic_mean_pptdorm$Sex==s])^2 * (s-1) +
                    (panic_coef$btempgrow2sex_p[post_draws[p]]) * (climgrow$tempgrow)^2 * (s-1) +
                    (panic_coef$btempdorm2sex_p[post_draws[p]]) * (panic_mean_tempdorm$tempdorm[panic_mean_tempdorm$Sex==s])^2 * (s-1)
                  
    )
    surv_sex_diff_post_grow[i,,p] <-  fem_s-male_s
    grow_sex_diff_post_grow[i,,p] <-  fem_g-male_g
    flow_sex_diff_post_grow[i,,p] <-  fem_f-male_f
    panic_sex_diff_post_grow[i,,p] <-  fem_p-male_p
  }
}


# define quantiles of posterior samples
sex_diff_surv_mean_grow <- sex_diff_grow_mean_grow <- sex_diff_flow_mean_grow <- sex_diff_panic_mean_grow <- matrix(NA,size_bin_num,nrow(climgrow))
  
  for(s in 1:size_bin_num){
    for(l in 1:nrow(climgrow)){
      sex_diff_surv_mean_grow[s,l] <- mean(surv_sex_diff_post_grow[s,l,],na.rm=TRUE)
      sex_diff_grow_mean_grow[s,l] <- mean(grow_sex_diff_post_grow[s,l,],na.rm=TRUE)
      sex_diff_flow_mean_grow[s,l] <- mean(flow_sex_diff_post_grow[s,l,],na.rm=TRUE)
      sex_diff_panic_mean_grow[s,l] <- mean(panic_sex_diff_post_grow[s,l,],na.rm=TRUE)
      
    }
  }

scal_breaks_diff_surv_growing <- c(seq(-0.3, 0.2, length.out = 101))
scal_breaks_diff_grow_growing <- c(seq(-30, 40, length.out = 101))
scal_breaks_diff_flow_growing <- c(seq(-0.2, 0.21, length.out = 101))
scal_breaks_diff_panic_growing <- c(seq(-3.6, 3, length.out = 101))

# library(wesanderson)
# pal <- wes_palette("Zissou1", 100, type = "continuous")
pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Vital_rate_growing_3D.pdf",width=9,height=8,useDingbats = F)
par(mar=c(5,5,2,3),mfrow=c(2,2))

mtrx_surv_growing <- matrix(as.vector(sex_diff_surv_mean_grow), nrow = 30, dimnames = list(pptgrow_seq,tempgrow_seq))
fields::image.plot(pptgrow_seq*sd(poar_2015_2016$pptgrow)+mean(poar_2015_2016$pptgrow),tempgrow_seq*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow),mtrx_surv_growing,col=cm.colors(100),xlab="Growing season precip",ylab="Growing season temp",main="",breaks = scal_breaks_diff_surv_growing,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,legend.mar = 4,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste(Delta," (F-M)")), side=3, font=3, line=0.3, cex=0.6)) 
contour(pptgrow_seq*sd(poar_2015_2016$pptgrow)+mean(poar_2015_2016$pptgrow),tempgrow_seq*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow),mtrx_surv_growing,add=T,labcex=0.75,col="black") 
mtext("Survival",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "A",side = 3, adj = 0,cex=1.2)

mtrx_grow_growing <- matrix(as.vector(sex_diff_grow_mean_grow), nrow = 30, dimnames = list(pptgrow_seq,tempgrow_seq))
fields::image.plot(pptgrow_seq*sd(poar_2015_2016$pptgrow)+mean(poar_2015_2016$pptgrow),tempgrow_seq*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow),mtrx_grow_growing,col=cm.colors(100),xlab="Growing season precip",ylab="Growing season temp",main="",breaks = scal_breaks_diff_grow_growing,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste(Delta," (F-M)")), side=3, font=2, line=0.3, cex=0.6)) 
contour(pptgrow_seq*sd(poar_2015_2016$pptgrow)+mean(poar_2015_2016$pptgrow),tempgrow_seq*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow),mtrx_grow_growing,add=T,labcex=0.75,col="black") 
mtext("Growth",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "B",side = 3, adj = 0,cex=1.2)

mtrx_flow_growing <- matrix(as.vector(sex_diff_flow_mean_grow), nrow = 30, dimnames = list(pptgrow_seq,tempgrow_seq))
fields::image.plot(pptgrow_seq*sd(poar_2015_2016$pptgrow)+mean(poar_2015_2016$pptgrow),tempgrow_seq*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow),mtrx_flow_growing,col=cm.colors(100),xlab="Growing season precip",ylab="Growing season temp",main="",breaks = scal_breaks_diff_flow_growing,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste(Delta," (F-M)")), side=3, font=2, line=0.3, cex=0.6)) 
contour(pptgrow_seq*sd(poar_2015_2016$pptgrow)+mean(poar_2015_2016$pptgrow),tempgrow_seq*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow),mtrx_flow_growing,add=T,labcex=0.75,col="black") 
mtext("Flowering",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "C",side = 3, adj = 0,cex=1.2)

mtrx_panic_growing <- matrix(as.vector(sex_diff_panic_mean_grow), nrow = 30, dimnames = list(pptgrow_seq,tempgrow_seq))
fields::image.plot(pptgrow_seq*sd(poar_2015_2016$pptgrow)+mean(poar_2015_2016$pptgrow),tempgrow_seq*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow),mtrx_panic_growing,col=cm.colors(100),xlab="Growing season temp",ylab="Growing season temp",main="",breaks = scal_breaks_diff_panic_growing,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste(Delta," (F-M)")), side=3, font=2, line=0.3, cex=0.6)) 
contour(pptgrow_seq*sd(poar_2015_2016$pptgrow)+mean(poar_2015_2016$pptgrow),tempgrow_seq*sd(poar_2015_2016$tempgrow)+mean(poar_2015_2016$tempgrow),mtrx_panic_growing,add=T,labcex=0.75,col="black") 
mtext("Panicles",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "D",side = 3, adj = 0,cex=1.2)

dev.off()

## Dormant season----
pptdorm_seq <- seq(min(poar_surv_binned$pptdorm),max(poar_surv_binned$pptdorm),length.out=30)
tempdorm_seq <- seq(min(poar_surv_binned$tempdorm),max(poar_surv_binned$tempdorm),length.out=30)
climdorm<-expand.grid(pptdorm=pptdorm_seq,tempdorm=tempdorm_seq)

n_post_draws <- 2000
post_draws <- sample.int(length(surv_coef$b0_s), n_post_draws)
surv_sex_diff_post_dorm <- grow_sex_diff_post_dorm <- flow_sex_diff_post_dorm <- panic_sex_diff_post_dorm <- array(NA,dim=c(size_bin_num,nrow(climdorm),n_post_draws))
for(p in 1:n_post_draws){
  for(i in 1:size_bin_num){
    s=1;      
    fem_s <- invlogit((surv_coef$b0_s[post_draws[p]]) + 
                        (surv_coef$bsize_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i]+
                        (surv_coef$bsizesex_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                        (surv_coef$bsex_s[post_draws[p]]) * (s-1) +
                        (surv_coef$bpptgrow_s[post_draws[p]]) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] +
                        (surv_coef$bpptdorm_s[post_draws[p]]) *  climdorm$pptdorm +
                        (surv_coef$btempgrow_s[post_draws[p]]) * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] +
                        (surv_coef$btempdorm_s[post_draws[p]]) * climdorm$tempdorm +
                        (surv_coef$bpptgrowsex_s[post_draws[p]]) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * (s-1) +
                        (surv_coef$bpptdormsex_s[post_draws[p]]) * climdorm$pptdorm * (s-1) +
                        (surv_coef$btempgrowsex_s[post_draws[p]]) * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] * (s-1) +
                        (surv_coef$btempdormsex_s[post_draws[p]]) * climdorm$tempdorm * (s-1) +
                        (surv_coef$btempdormpptdorm_s[post_draws[p]]) * climdorm$pptdorm * climdorm$tempdorm +
                        (surv_coef$btempgrowpptgrow_s[post_draws[p]]) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] +
                        (surv_coef$btempdormpptdormsex_s[post_draws[p]]) * climdorm$tempdorm * climdorm$pptdorm * (s-1) +
                        (surv_coef$btempgrowpptgrowsex_s[post_draws[p]]) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] * (s-1) +
                        (surv_coef$bpptgrow2_s[post_draws[p]]) * (surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s])^2 +
                        (surv_coef$bpptdorm2_s[post_draws[p]]) * (climdorm$pptdorm)^2 +
                        (surv_coef$btempgrow2_s[post_draws[p]]) * (surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s])^2 +
                        (surv_coef$btempdorm2_s[post_draws[p]]) * (climdorm$tempdorm)^2 +
                        (surv_coef$bpptgrow2sex_s[post_draws[p]]) * (surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s])^2 * (s-1) +
                        (surv_coef$bpptdorm2sex_s[post_draws[p]]) * (climdorm$pptdorm)^2 * (s-1) +
                        (surv_coef$btempgrow2sex_s[post_draws[p]]) * (surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s])^2 * (s-1) +
                        (surv_coef$btempdorm2sex_s[post_draws[p]]) * (climdorm$tempdorm)^2 * (s-1)
    )
    
    fem_g <- exp((grow_coef$b0_g[post_draws[p]]) + 
                   (grow_coef$bsize_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i]+
                   (grow_coef$bsizesex_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                   (grow_coef$bsex_g[post_draws[p]]) * (s-1) +
                   (grow_coef$bpptgrow_g[post_draws[p]]) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] +
                   (grow_coef$bpptdorm_g[post_draws[p]]) *  climdorm$pptdorm +
                   (grow_coef$btempgrow_g[post_draws[p]]) * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] +
                   (grow_coef$btempdorm_g[post_draws[p]]) * climdorm$tempdorm +
                   (grow_coef$bpptgrowsex_g[post_draws[p]]) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * (s-1) +
                   (grow_coef$bpptdormsex_g[post_draws[p]]) * climdorm$pptdorm * (s-1) +
                   (grow_coef$btempgrowsex_g[post_draws[p]]) * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] * (s-1) +
                   (grow_coef$btempdormsex_g[post_draws[p]]) * climdorm$tempdorm * (s-1) +
                   (grow_coef$btempdormpptdorm_g[post_draws[p]]) * climdorm$pptdorm * climdorm$tempdorm +
                   (grow_coef$btempgrowpptgrow_g[post_draws[p]]) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] +
                   (grow_coef$btempdormpptdormsex_g[post_draws[p]]) * climdorm$tempdorm * climdorm$pptdorm * (s-1) +
                   (grow_coef$btempgrowpptgrowsex_g[post_draws[p]]) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] * (s-1) +
                   (grow_coef$bpptgrow2_g[post_draws[p]]) * (grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s])^2 +
                   (grow_coef$bpptdorm2_g[post_draws[p]]) * (climdorm$pptdorm)^2 +
                   (grow_coef$btempgrow2_g[post_draws[p]]) * (grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s])^2 +
                   (grow_coef$btempdorm2_g[post_draws[p]]) * (climdorm$tempdorm)^2 +
                   (grow_coef$bpptgrow2sex_g[post_draws[p]]) * (grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s])^2 * (s-1) +
                   (grow_coef$bpptdorm2sex_g[post_draws[p]]) * (climdorm$pptdorm)^2 * (s-1) +
                   (grow_coef$btempgrow2sex_g[post_draws[p]]) * (grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s])^2 * (s-1) +
                   (grow_coef$btempdorm2sex_g[post_draws[p]]) * (climdorm$tempdorm)^2 * (s-1)
    )
    
    fem_f <- invlogit((flow_coef$b0_f[post_draws[p]]) + 
                        (flow_coef$bsize_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i]+
                        (flow_coef$bsizesex_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                        (flow_coef$bsex_f[post_draws[p]]) * (s-1) +
                        (flow_coef$bpptgrow_f[post_draws[p]]) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] +
                        (flow_coef$bpptdorm_f[post_draws[p]]) *  climdorm$pptdorm +
                        (flow_coef$btempgrow_f[post_draws[p]]) * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] +
                        (flow_coef$btempdorm_f[post_draws[p]]) * climdorm$tempdorm +
                        (flow_coef$bpptgrowsex_f[post_draws[p]]) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * (s-1) +
                        (flow_coef$bpptdormsex_f[post_draws[p]]) * climdorm$pptdorm * (s-1) +
                        (flow_coef$btempgrowsex_f[post_draws[p]]) * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] * (s-1) +
                        (flow_coef$btempdormsex_f[post_draws[p]]) * climdorm$tempdorm * (s-1) +
                        (flow_coef$btempdormpptdorm_f[post_draws[p]]) * climdorm$pptdorm * climdorm$tempdorm +
                        (flow_coef$btempgrowpptgrow_f[post_draws[p]]) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] +
                        (flow_coef$btempdormpptdormsex_f[post_draws[p]]) * climdorm$tempdorm * climdorm$pptdorm  * (s-1) +
                        (flow_coef$btempgrowpptgrowsex_f[post_draws[p]]) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] * (s-1) +
                        (flow_coef$bpptgrow2_f[post_draws[p]]) * (flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s])^2 +
                        (flow_coef$bpptdorm2_f[post_draws[p]]) * (climdorm$pptdorm)^2 +
                        (flow_coef$btempgrow2_f[post_draws[p]]) * (flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s])^2 +
                        (flow_coef$btempdorm2_f[post_draws[p]]) * (climdorm$tempdorm)^2 +
                        (flow_coef$bpptgrow2sex_f[post_draws[p]]) * (flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s])^2 * (s-1) +
                        (flow_coef$bpptdorm2sex_f[post_draws[p]]) * (climdorm$pptdorm)^2 * (s-1) +
                        (flow_coef$btempgrow2sex_f[post_draws[p]]) * (flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s])^2 * (s-1) +
                        (flow_coef$btempdorm2sex_f[post_draws[p]]) * (climdorm$tempdorm)^2 * (s-1)
    )
    
    fem_p <- exp((panic_coef$b0_p[post_draws[p]]) + 
                   (panic_coef$bsize_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i]+
                   (panic_coef$bsizesex_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                   (panic_coef$bsex_p[post_draws[p]]) * (s-1) +
                   (panic_coef$bpptgrow_p[post_draws[p]]) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] +
                   (panic_coef$bpptdorm_p[post_draws[p]]) *  climdorm$pptdorm +
                   (panic_coef$btempgrow_p[post_draws[p]]) * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] +
                   (panic_coef$btempdorm_p[post_draws[p]]) * climdorm$tempdorm +
                   (panic_coef$bpptgrowsex_p[post_draws[p]]) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] * (s-1) +
                   (panic_coef$bpptdormsex_p[post_draws[p]]) * climdorm$pptdorm * (s-1) +
                   (panic_coef$btempgrowsex_p[post_draws[p]]) * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] * (s-1) +
                   (panic_coef$btempdormsex_p[post_draws[p]]) * climdorm$tempdorm * (s-1) +
                   (panic_coef$btempdormpptdorm_p[post_draws[p]]) * climdorm$pptdorm * climdorm$tempdorm +
                   (panic_coef$btempgrowpptgrow_p[post_draws[p]]) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] +
                   (panic_coef$btempdormpptdormsex_p[post_draws[p]]) * climdorm$tempdorm * climdorm$pptdorm * (s-1) +
                   (panic_coef$btempgrowpptgrowsex_p[post_draws[p]]) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] * (s-1) +
                   (panic_coef$bpptgrow2_p[post_draws[p]]) * (panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s])^2 +
                   (panic_coef$bpptdorm2_p[post_draws[p]]) * (climdorm$pptdorm)^2 +
                   (panic_coef$btempgrow2_p[post_draws[p]]) * (panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s])^2 +
                   (panic_coef$btempdorm2_p[post_draws[p]]) * (climdorm$tempdorm)^2 +
                   (panic_coef$bpptgrow2sex_p[post_draws[p]]) * (panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s])^2 * (s-1) +
                   (panic_coef$bpptdorm2sex_p[post_draws[p]]) * (climdorm$pptdorm)^2 * (s-1) +
                   (panic_coef$btempgrow2sex_p[post_draws[p]]) * (panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s])^2 * (s-1) +
                   (panic_coef$btempdorm2sex_p[post_draws[p]]) * (climdorm$tempdorm)^2 * (s-1)
    )
    
    s=2;       
    male_s <- invlogit((surv_coef$b0_s[post_draws[p]]) + 
                         (surv_coef$bsize_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i]+
                         (surv_coef$bsizesex_s[post_draws[p]]) * surv_mean_sizes$size[surv_mean_sizes$Sex==s & surv_mean_sizes$size_bin==i] * (s-1) +
                         (surv_coef$bsex_s[post_draws[p]]) * (s-1) +
                         (surv_coef$bpptgrow_s[post_draws[p]]) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] +
                         (surv_coef$bpptdorm_s[post_draws[p]]) *  climdorm$pptdorm +
                         (surv_coef$btempgrow_s[post_draws[p]]) * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] +
                         (surv_coef$btempdorm_s[post_draws[p]]) * climdorm$tempdorm +
                         (surv_coef$bpptgrowsex_s[post_draws[p]]) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * (s-1) +
                         (surv_coef$bpptdormsex_s[post_draws[p]]) * climdorm$pptdorm * (s-1) +
                         (surv_coef$btempgrowsex_s[post_draws[p]]) * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] * (s-1) +
                         (surv_coef$btempdormsex_s[post_draws[p]]) * climdorm$tempdorm * (s-1) +
                         (surv_coef$btempdormpptdorm_s[post_draws[p]]) * climdorm$pptdorm * climdorm$tempdorm +
                         (surv_coef$btempgrowpptgrow_s[post_draws[p]]) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] +
                         (surv_coef$btempdormpptdormsex_s[post_draws[p]]) * climdorm$tempdorm * climdorm$pptdorm * (s-1) +
                         (surv_coef$btempgrowpptgrowsex_s[post_draws[p]]) * surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s] * surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s] * (s-1) +
                         (surv_coef$bpptgrow2_s[post_draws[p]]) * (surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s])^2 +
                         (surv_coef$bpptdorm2_s[post_draws[p]]) * (climdorm$pptdorm)^2 +
                         (surv_coef$btempgrow2_s[post_draws[p]]) * (surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s])^2 +
                         (surv_coef$btempdorm2_s[post_draws[p]]) * (climdorm$tempdorm)^2 +
                         (surv_coef$bpptgrow2sex_s[post_draws[p]]) * (surv_mean_pptgrow$pptgrow[surv_mean_pptgrow$Sex==s])^2 * (s-1) +
                         (surv_coef$bpptdorm2sex_s[post_draws[p]]) * (climdorm$pptdorm)^2 * (s-1) +
                         (surv_coef$btempgrow2sex_s[post_draws[p]]) * (surv_mean_tempgrow$tempgrow[surv_mean_tempgrow$Sex==s])^2 * (s-1) +
                         (surv_coef$btempdorm2sex_s[post_draws[p]]) * (climdorm$tempdorm)^2 * (s-1)
    )
    
    male_g <- exp((grow_coef$b0_g[post_draws[p]]) + 
                    (grow_coef$bsize_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i]+
                    (grow_coef$bsizesex_g[post_draws[p]]) * grow_mean_sizes$size[grow_mean_sizes$Sex==s & grow_mean_sizes$size_bin==i] * (s-1) +
                    (grow_coef$bsex_g[post_draws[p]]) * (s-1) +
                    (grow_coef$bpptgrow_g[post_draws[p]]) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] +
                    (grow_coef$bpptdorm_g[post_draws[p]]) *  climdorm$pptdorm +
                    (grow_coef$btempgrow_g[post_draws[p]]) * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] +
                    (grow_coef$btempdorm_g[post_draws[p]]) * climdorm$tempdorm +
                    (grow_coef$bpptgrowsex_g[post_draws[p]]) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * (s-1) +
                    (grow_coef$bpptdormsex_g[post_draws[p]]) * climdorm$pptdorm * (s-1) +
                    (grow_coef$btempgrowsex_g[post_draws[p]]) * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] * (s-1) +
                    (grow_coef$btempdormsex_g[post_draws[p]]) * climdorm$tempdorm * (s-1) +
                    (grow_coef$btempdormpptdorm_g[post_draws[p]]) * climdorm$pptdorm * climdorm$tempdorm +
                    (grow_coef$btempgrowpptgrow_g[post_draws[p]]) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] +
                    (grow_coef$btempdormpptdormsex_g[post_draws[p]]) * climdorm$tempdorm * climdorm$pptdorm * (s-1) +
                    (grow_coef$btempgrowpptgrowsex_g[post_draws[p]]) * grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s] * grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s] * (s-1) +
                    (grow_coef$bpptgrow2_g[post_draws[p]]) * (grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s])^2 +
                    (grow_coef$bpptdorm2_g[post_draws[p]]) * (climdorm$pptdorm)^2 +
                    (grow_coef$btempgrow2_g[post_draws[p]]) * (grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s])^2 +
                    (grow_coef$btempdorm2_g[post_draws[p]]) * (climdorm$tempdorm)^2 +
                    (grow_coef$bpptgrow2sex_g[post_draws[p]]) * (grow_mean_pptgrow$pptgrow[grow_mean_pptgrow$Sex==s])^2 * (s-1) +
                    (grow_coef$bpptdorm2sex_g[post_draws[p]]) * (climdorm$pptdorm)^2 * (s-1) +
                    (grow_coef$btempgrow2sex_g[post_draws[p]]) * (grow_mean_tempgrow$tempgrow[grow_mean_tempgrow$Sex==s])^2 * (s-1) +
                    (grow_coef$btempdorm2sex_g[post_draws[p]]) * (climdorm$tempdorm)^2 * (s-1)
    )
    
    male_f <- invlogit((flow_coef$b0_f[post_draws[p]]) + 
                         (flow_coef$bsize_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i]+
                         (flow_coef$bsizesex_f[post_draws[p]]) * flow_mean_sizes$size[flow_mean_sizes$Sex==s & flow_mean_sizes$size_bin==i] * (s-1) +
                         (flow_coef$bsex_f[post_draws[p]]) * (s-1) +
                         (flow_coef$bpptgrow_f[post_draws[p]]) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] +
                         (flow_coef$bpptdorm_f[post_draws[p]]) *  climdorm$pptdorm +
                         (flow_coef$btempgrow_f[post_draws[p]]) * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] +
                         (flow_coef$btempdorm_f[post_draws[p]]) * climdorm$tempdorm +
                         (flow_coef$bpptgrowsex_f[post_draws[p]]) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * (s-1) +
                         (flow_coef$bpptdormsex_f[post_draws[p]]) * climdorm$pptdorm * (s-1) +
                         (flow_coef$btempgrowsex_f[post_draws[p]]) * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] * (s-1) +
                         (flow_coef$btempdormsex_f[post_draws[p]]) * climdorm$tempdorm * (s-1) +
                         (flow_coef$btempdormpptdorm_f[post_draws[p]]) * climdorm$pptdorm * climdorm$tempdorm +
                         (flow_coef$btempgrowpptgrow_f[post_draws[p]]) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] +
                         (flow_coef$btempdormpptdormsex_f[post_draws[p]]) * climdorm$tempdorm * climdorm$pptdorm  * (s-1) +
                         (flow_coef$btempgrowpptgrowsex_f[post_draws[p]]) * flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s] * flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s] * (s-1) +
                         (flow_coef$bpptgrow2_f[post_draws[p]]) * (flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s])^2 +
                         (flow_coef$bpptdorm2_f[post_draws[p]]) * (climdorm$pptdorm)^2 +
                         (flow_coef$btempgrow2_f[post_draws[p]]) * (flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s])^2 +
                         (flow_coef$btempdorm2_f[post_draws[p]]) * (climdorm$tempdorm)^2 +
                         (flow_coef$bpptgrow2sex_f[post_draws[p]]) * (flow_mean_pptgrow$pptgrow[flow_mean_pptgrow$Sex==s])^2 * (s-1) +
                         (flow_coef$bpptdorm2sex_f[post_draws[p]]) * (climdorm$pptdorm)^2 * (s-1) +
                         (flow_coef$btempgrow2sex_f[post_draws[p]]) * (flow_mean_tempgrow$tempgrow[flow_mean_tempgrow$Sex==s])^2 * (s-1) +
                         (flow_coef$btempdorm2sex_f[post_draws[p]]) * (climdorm$tempdorm)^2 * (s-1)
    )
    
    
    male_p <- exp((panic_coef$b0_p[post_draws[p]]) + 
                    (panic_coef$bsize_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i]+
                    (panic_coef$bsizesex_p[post_draws[p]]) * panic_mean_sizes$size[panic_mean_sizes$Sex==s & panic_mean_sizes$size_bin==i] * (s-1) +
                    (panic_coef$bsex_p[post_draws[p]]) * (s-1) +
                    (panic_coef$bpptgrow_p[post_draws[p]]) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] +
                    (panic_coef$bpptdorm_p[post_draws[p]]) *  climdorm$pptdorm +
                    (panic_coef$btempgrow_p[post_draws[p]]) * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] +
                    (panic_coef$btempdorm_p[post_draws[p]]) * climdorm$tempdorm +
                    (panic_coef$bpptgrowsex_p[post_draws[p]]) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] * (s-1) +
                    (panic_coef$bpptdormsex_p[post_draws[p]]) * climdorm$pptdorm * (s-1) +
                    (panic_coef$btempgrowsex_p[post_draws[p]]) * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] * (s-1) +
                    (panic_coef$btempdormsex_p[post_draws[p]]) * climdorm$tempdorm * (s-1) +
                    (panic_coef$btempdormpptdorm_p[post_draws[p]]) * climdorm$pptdorm * climdorm$tempdorm +
                    (panic_coef$btempgrowpptgrow_p[post_draws[p]]) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] +
                    (panic_coef$btempdormpptdormsex_p[post_draws[p]]) * climdorm$tempdorm * climdorm$pptdorm * (s-1) +
                    (panic_coef$btempgrowpptgrowsex_p[post_draws[p]]) * panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s] * panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s] * (s-1) +
                    (panic_coef$bpptgrow2_p[post_draws[p]]) * (panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s])^2 +
                    (panic_coef$bpptdorm2_p[post_draws[p]]) * (climdorm$pptdorm)^2 +
                    (panic_coef$btempgrow2_p[post_draws[p]]) * (panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s])^2 +
                    (panic_coef$btempdorm2_p[post_draws[p]]) * (climdorm$pptdorm)^2 +
                    (panic_coef$bpptgrow2sex_p[post_draws[p]]) * (panic_mean_pptgrow$pptgrow[panic_mean_pptgrow$Sex==s])^2 * (s-1) +
                    (panic_coef$bpptdorm2sex_p[post_draws[p]]) * (climdorm$pptdorm)^2 * (s-1) +
                    (panic_coef$btempgrow2sex_p[post_draws[p]]) * (panic_mean_tempgrow$tempgrow[panic_mean_tempgrow$Sex==s])^2 * (s-1) +
                    (panic_coef$btempdorm2sex_p[post_draws[p]]) * (climdorm$pptdorm)^2 * (s-1)
    )
    surv_sex_diff_post_dorm[i,,p] <-  fem_s-male_s
    grow_sex_diff_post_dorm[i,,p] <-  fem_g-male_g
    flow_sex_diff_post_dorm[i,,p] <-  fem_f-male_f
    panic_sex_diff_post_dorm[i,,p] <-  fem_p-male_p
  }
}

#define quantiles of posterior samples
sex_diff_surv_mean_dorm <- sex_diff_grow_mean_dorm <- sex_diff_flow_mean_dorm <- sex_diff_panic_mean_dorm <- matrix(NA,size_bin_num,nrow(climdorm))
for(s in 1:size_bin_num){
  for(l in 1:nrow(climdorm)){
    sex_diff_surv_mean_dorm[s,l] <- mean(surv_sex_diff_post_dorm[s,l,],na.rm=TRUE)
    sex_diff_grow_mean_dorm[s,l] <- mean(grow_sex_diff_post_dorm[s,l,],na.rm=TRUE)
    sex_diff_flow_mean_dorm[s,l] <- mean(flow_sex_diff_post_dorm[s,l,],na.rm=TRUE)
    sex_diff_panic_mean_dorm[s,l] <- mean(panic_sex_diff_post_dorm[s,l,],na.rm=TRUE)
    
  }
}

scal_breaks_diff_surv_dormant <- c(seq(-0.6, 0.6, length.out = 101))
scal_breaks_diff_grow_dormant <- c(seq(-85,87, length.out = 101))
scal_breaks_diff_flow_dormant <- c(seq(-0.2, 0.25, length.out = 101))
scal_breaks_diff_panic_dormant <- c(seq(-3.8, 3, length.out = 101))


pdf("/Users/jm200/Library/CloudStorage/Dropbox/Miller Lab/github/POAR-Forecasting/Manuscript/Figures/Vital_rate_dormant_3D.pdf",width=9,height=8,useDingbats = F)
par(mar=c(5,5,2,3),mfrow=c(2,2))
mtrx_surv_dormant <- matrix(as.vector(sex_diff_surv_mean_dorm), nrow = 30, dimnames = list(pptdorm_seq,tempdorm_seq))
fields::image.plot(pptdorm_seq*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm),tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mtrx_surv_dormant,col=cm.colors(100),xlab="Dormant season precip",ylab="Dormant season temp",main="",breaks = scal_breaks_diff_surv_dormant,cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,legend.mar = 4,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste(Delta," (F-M)")), side=3, font=3, line=0.3, cex=0.6))
contour(pptdorm_seq*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm),tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mtrx_surv_dormant,add=T,labcex=0.75,col="black") 
mtext("Survival",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "A",side = 3, adj = 0,cex=1.2)
mtrx_grow_dormant <- matrix(as.vector(sex_diff_grow_mean_dorm), nrow = 30, dimnames = list(pptdorm_seq,tempdorm_seq))
fields::image.plot(pptdorm_seq*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm),tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mtrx_grow_dormant,col=cm.colors(100),xlab="Dorm. season precip",ylab="Dormant season temp",main="",breaks = scal_breaks_diff_grow_dormant, cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,legend.mar = 4,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste(Delta," (F-M)")), side=3, font=3, line=0.3, cex=0.6)) 
contour(pptdorm_seq*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm),tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mtrx_grow_dormant,add=T,labcex=0.75,col="black") 
mtext("Growth",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "B",side = 3, adj = 0,cex=1.2)
mtrx_flow_dormant <- matrix(as.vector(sex_diff_flow_mean_dorm), nrow = 30, dimnames = list(pptdorm_seq,tempdorm_seq))
fields::image.plot(pptdorm_seq*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm),tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mtrx_flow_dormant,col=cm.colors(100),xlab="Dormant season precip",ylab="Dormant season temp",main="",breaks = scal_breaks_diff_flow_dormant, cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,legend.mar = 4,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste(Delta," (F-M)")), side=3, font=3, line=0.3, cex=0.6)) 
contour(pptdorm_seq*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm),tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mtrx_flow_dormant,add=T,labcex=0.75,col="black") 
mtext("Flowering",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "C",side = 3, adj = 0,cex=1.2)
mtrx_panic_dormant <- matrix(as.vector(sex_diff_panic_mean_dorm), nrow = 30, dimnames = list(pptdorm_seq,tempdorm_seq))
fields::image.plot(pptdorm_seq*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm),tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mtrx_panic_dormant,col=cm.colors(100),xlab="Dormant season precip",ylab="Dormant season temp",main="",breaks = scal_breaks_diff_panic_dormant, cex.lab=1.2,
                   legend.width=1, legend.shrink=0.75,legend.mar = 4,
                   axis.args=list(cex.axis=0.6),
                   legend.args=list(text=expression(paste(Delta," (F-M)")), side=3, font=3, line=0.3, cex=0.6)) 
contour(pptdorm_seq*sd(poar_2015_2016$pptdorm)+mean(poar_2015_2016$pptdorm),tempdorm_seq*sd(poar_2015_2016$tempdorm) + mean(poar_2015_2016$tempdorm),mtrx_panic_dormant,add=T,labcex=0.75,col="black") 
mtext("Panicles",side = 3, adj = 0.5,cex=1.2,line=0.3)
mtext( "D",side = 3, adj = 0,cex=1.2)
dev.off()

