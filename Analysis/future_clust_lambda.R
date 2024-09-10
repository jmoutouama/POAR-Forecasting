# The goal of this code is to facilitate the estimation of population growth rate estimation across Texas, Oklahoma and Kansas using the Rice Super computing (NOTS)
library(tidyverse)
library(actuar)
set.seed(13)

# Benchmark the time
init_t  <- Sys.time()

# arguments from command line
args    <- commandArgs(trailingOnly = TRUE)

# tests
print('let us see')
print(args)
class(args)
class(args[1])
print( getwd() )

# identify the site
job_id      <- as.numeric(args[1])
future_prj  <- as.character( args[2] )

print( job_id )

# design model
design_df <- expand.grid( clim_id = 1:3972,
                          post_id = 1:300 )

l             <- design_df[job_id,]$clim_id
p             <- design_df[job_id,]$post_id
max_yrs       <- 30

# quote a series of bare names
quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}


# read data --------------------------------------------------------------------

# seedling data
sdlg_surv               <- read.csv( paste0( '/gpfs1/data/idiv_knight/sApropos/',
                                             'poar_climate/sdlg_surv.csv' ) 
                                     )

# seedling data
poar.clim_seasonal      <- read.csv( paste0( '/gpfs1/data/idiv_knight/sApropos/',
                                             'poar_climate/poar.clim_seasonal.csv' ) 
                                     )

# current climate
future_file             <- paste0('climate_', future_prj, '_values.csv') 
climate_values          <- read.csv( paste0( '/gpfs1/data/idiv_knight/sApropos/',
                                             'poar_climate/', future_file ) 
                                     )

# full Bayeisan model
post_draws              <- read.csv( paste0( '/gpfs1/data/idiv_knight/sApropos/',
                                             'poar_climate/post_draws.csv' )
                                     ) %>% .[,1]

# full Bayeisan model
fit_full                <- readRDS( paste0( '/gpfs1/data/idiv_knight/sApropos/',
                                            'poar_climate/full_model.RDS' ) )


# vital rate posteriors --------------------------------------------------------

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

# load MPM functions

# vital rate and megamatrix functions ---------------------------------------------------
invlogit<-function(x){exp(x)/(1+exp(x))}

#SURVIVAL AT SIZE X.
sx<-function(x,params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,surv_perturb=0){
  surv_mean<-params$surv_mu + 
    params$surv_size*log(x) + 
    params$surv_pptgrow*zpptgrow + 
    params$surv_pptdorm*zpptdorm+ 
    params$surv_tempgrow*ztempgrow + 
    params$surv_tempdorm*ztempdorm + 
    params$surv_tempdorm_pptdorm*ztempdorm*zpptdorm + 
    params$surv_tempgrow_pptgrow*ztempgrow*zpptgrow + 
    params$surv_pptgrow2*(zpptgrow^2) + 
    params$surv_pptdorm2*(zpptdorm^2) + 
    params$surv_tempgrow2*(ztempgrow^2) + 
    params$surv_tempdorm2*(ztempdorm^2) +
    rfx["surv","site"] + rfx["surv","block"] + rfx["surv","source"]
  return(invlogit(surv_mean)+ surv_perturb)
}


#PROBABILITY OF GROWTH FROM SIZE X TO Y
#This function truncates the density asscociation with x==0 and x>x.max
gxy<-function(x,y,params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,grow_perturb=0){
  grow_mean<-exp(params$grow_mu + 
                   params$grow_size*log(x)+ 
                   params$grow_pptgrow*zpptgrow + 
                   params$grow_pptdorm*zpptdorm + 
                   params$grow_tempgrow*ztempgrow + 
                   params$grow_tempdorm*ztempdorm + 
                   params$grow_tempdorm_pptdorm*ztempdorm*zpptdorm + 
                   params$grow_tempgrow_pptgrow*ztempgrow*zpptgrow + 
                   params$grow_pptgrow2*(zpptgrow^2) + 
                   params$grow_pptdorm2*(zpptdorm^2) + 
                   params$grow_tempgrow2*(ztempgrow^2) + 
                   params$grow_tempdorm2*(ztempdorm^2) +
                   rfx["grow","site"] + rfx["grow","block"] + rfx["grow","source"]) + grow_perturb
  grow<-actuar::dpoisinvgauss(x=y,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  grow<-ifelse(is.nan(grow) | is.infinite(grow),0,grow)
  truncLower<-actuar::dpoisinvgauss(x=0,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  truncUpper<-actuar::dpoisinvgauss(x=params$max_size:10000,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  truncUpper<-sum(ifelse(is.nan(truncUpper) | is.infinite(truncUpper),0,truncUpper))
  return(grow/(1-(truncLower+truncUpper)))
}

#SURVIVAL*GROWTH 
pxy<-function(x,y,params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,surv_perturb=0,grow_perturb=0){
  sx(x,params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,surv_perturb) * gxy(x,y,params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,grow_perturb)
}


# PROBABILITY OF FLOWERING
pfx<-function(x,params, zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,flow_perturb=0){
  flow_mean<-params$flow_mu + 
    params$flow_size*log(x) + 
    params$flow_pptgrow*zpptgrow + 
    params$flow_pptdorm*zpptdorm + 
    params$flow_tempgrow*ztempgrow + 
    params$flow_tempdorm*ztempdorm + 
    params$flow_tempdorm_pptdorm*ztempdorm*zpptdorm + 
    params$flow_tempgrow_pptgrow*ztempgrow*zpptgrow + 
    params$flow_pptgrow2*(zpptgrow^2) + 
    params$flow_pptdorm2*(zpptdorm^2) + 
    params$flow_tempgrow2*(ztempgrow^2) + 
    params$flow_tempdorm2*(ztempdorm^2) +
    rfx["flow","site"] + rfx["flow","block"] + rfx["flow","source"]
  return(invlogit(flow_mean) + flow_perturb )
}

#NUMBER OF PANICLES
nfx<-function(x,params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,fert_perturb=0){
  # x_scaled_p=(x-38.78495)/58.04579
  panic_mean<-params$panic_mu + 
    params$panic_size*log(x) + 
    params$panic_pptgrow*zpptgrow + 
    params$panic_pptdorm*zpptdorm + 
    params$panic_tempgrow*ztempgrow + 
    params$panic_tempdorm*ztempdorm + 
    params$panic_tempdorm_pptdorm*ztempdorm*zpptdorm + 
    params$panic_tempgrow_pptgrow*ztempgrow*zpptgrow + 
    params$panic_pptgrow2*(zpptgrow^2) + 
    params$panic_pptdorm2*(zpptdorm^2) + 
    params$panic_tempgrow2*(ztempgrow^2) + 
    params$panic_tempdorm2*(ztempdorm^2) +
    rfx["panic","site"] + rfx["panic","block"] + rfx["panic","source"]
  return(exp(panic_mean) + fert_perturb)
}


#SEED VIABILITY
viab<-function(params,twosex,OSR=NULL,viab_perturb=0){
  if(twosex==F){return(params$v0 + viab_perturb)}
  if(twosex==T){return((params$v0 * (1 - OSR ^ params$a_v))+ viab_perturb)}
}

#FERTILITY--returns number of recruits
##Female offspring
fertx_F<-function(x,params,rfx,zpptgrow,zpptdorm,ztempgrow,ztempdorm,twosex,OSR=NULL,flow_perturb=0,fert_perturb=0,viab_perturb=0){
  seedlings<-pfx(x,params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,flow_perturb)*nfx(x,params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,fert_perturb)*params$ov_per_inf*viab(params,twosex,OSR,viab_perturb)*params$germ*params$PSR
  return(seedlings)
}


##Male offspring
fertx_M<-function(x,params,rfx,zpptgrow,zpptdorm,ztempgrow,ztempdorm,twosex,OSR=NULL,flow_perturb=0,fert_perturb=0,viab_perturb=0){
  seedlings<-pfx(x,params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,flow_perturb)*nfx(x,params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,fert_perturb)*params$ov_per_inf*viab(params,twosex,OSR,viab_perturb)*params$germ*(1-params$PSR)
  return(seedlings)
}


megamatrix_delay<-function(F_params,M_params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,twosex,OSR=NULL,rfx,
                           surv_perturb=0,grow_perturb=0,flow_perturb=0,fert_perturb=0,viab_perturb=0){  
  matdim<-F_params$max_size       
  y<-1:F_params$max_size
  
  ## F-to-F (growth/survival transition)
  F.Tmat<-matrix(0,matdim+1,matdim+1)
  F.Tmat[2:(matdim+1),2:(matdim+1)]<-t(outer(y,y,pxy,params=F_params,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,rfx=rfx,surv_perturb=surv_perturb, grow_perturb=grow_perturb))
  # surviving seedlings emerge in continuous population
  F.Tmat[2:(matdim+1),1] <- gxy(x=1,y=y,params=F_params,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,rfx=rfx,grow_perturb=grow_perturb) * (M_params$sdlg_surv + surv_perturb)
  
  # F-to-F Fertility transition
  F.Fmat<-matrix(0,matdim+1,matdim+1)
  # seedlings in top row
  F.Fmat[1,2:(matdim+1)]<-fertx_F(x=y,params=F_params,rfx=rfx,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,twosex=twosex,OSR=OSR,flow_perturb=flow_perturb,fert_perturb=fert_perturb,viab_perturb=viab_perturb) 
  
  ## M-to-M (growth/survival transition)
  M.Tmat<-matrix(0,matdim+1,matdim+1)
  M.Tmat[2:(matdim+1),2:(matdim+1)]<-t(outer(y,y,pxy,params=M_params,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,rfx=rfx,surv_perturb=surv_perturb, grow_perturb=grow_perturb))
  M.Tmat[2:(matdim+1),1] <- gxy(x=1,y=y,params=M_params,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,rfx=rfx,grow_perturb=grow_perturb) * (M_params$sdlg_surv + surv_perturb)
  
  # F-to-M Fertility transition
  M.Fmat<-matrix(0,matdim+1,matdim+1)
  M.Fmat[1,2:(matdim+1)]<-fertx_M(x=y,params=F_params,rfx=rfx,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,twosex=twosex,OSR=OSR,flow_perturb=flow_perturb,fert_perturb=fert_perturb,viab_perturb=viab_perturb) 
  
  #M-to-F
  zero.mat<-matrix(0,matdim+1,matdim+1)
  
  # Put it all together as a megamatrix
  MEGAmat<-cbind(rbind(F.Tmat+F.Fmat,  ##Female growth/survival + recruitment[1,1]
                       M.Fmat), ##Male recruitment [2,1]
                 rbind(zero.mat,   ##Females from males [1,2]
                       M.Tmat))   ##Male growth/survival
  
  return(list(MEGAmat=MEGAmat,F.Tmat=F.Tmat,M.Tmat=M.Tmat,y=y))
}


# Analysis of 2sex model --------------------------------------------------

# this needs to be done by simulation
lambdaSim_delay<-function(F_params,M_params,zpptgrow,zpptdorm,ztempgrow,ztempdorm,rfx,max.yrs,
                          grow_perturb=0,surv_perturb=0,flow_perturb=0,fert_perturb=0,viab_perturb=0){
  matdim<-F_params$max_size       
  y<-1:F_params$max_size
  lambdatracker      <- rep(0,max.yrs)
  OSRtracker   <- rep(0,max.yrs)
  SRtracker   <- rep(0,max.yrs)
  n0            <- rep(1/((matdim+1)*2),((matdim+1)*2))
  
  for(t in 1:max.yrs){
    ##Estimate panicle SR
    flowering_females<-n0[2:(matdim+1)]*pfx(x=y,param=F_params,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,rfx=rfx,flow_perturb=flow_perturb) ## scalar multiplication to weight females by flowering prob
    F_panicles<-flowering_females%*%nfx(x=y,param=F_params,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,rfx=rfx,fert_perturb=fert_perturb) ##Vector operation to sum female panicles
    flowering_males<-n0[(matdim+3):((matdim+1)*2)]*pfx(x=y,param=M_params,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,rfx=rfx,flow_perturb=flow_perturb)
    M_panicles<-flowering_males %*% nfx(x=y,param=M_params,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,rfx=rfx,fert_perturb=fert_perturb)
    OSRtracker[t]<-F_panicles/(F_panicles+M_panicles) ##Panicle sex ratio (proportion female)
    SRtracker[t]<-sum(n0[1:(matdim+1)])
    #assmble matrix
    MEGAmat<-megamatrix_delay(F_params=F_params,M_params=M_params,zpptgrow=zpptgrow,zpptdorm=zpptdorm,ztempgrow=ztempgrow,ztempdorm=ztempdorm,twosex=T,OSR=OSRtracker[t],rfx=rfx,
                              grow_perturb=grow_perturb,surv_perturb=surv_perturb,flow_perturb=flow_perturb,fert_perturb=fert_perturb,viab_perturb=viab_perturb
    )$MEGAmat
    n0 <- MEGAmat[,] %*% n0
    N  <- sum(n0)
    lambdatracker[t]<-N
    n0 <-n0/N
  }
  return(list(lambdatracker=lambdatracker,SRtracker=SRtracker,OSRtracker=OSRtracker,n0=n0))
}


# RFX fun 
rfx_fun <- function(site_tau_s=0,block_tau_s=0,source_tau_s=0,
                    site_tau_g=0,block_tau_g=0,source_tau_g=0,
                    site_tau_f=0,block_tau_f=0,source_tau_f=0,
                    site_tau_p=0,block_tau_p=0,source_tau_p=0){
  rfx <- data.frame(site = rnorm(4,0,c(site_tau_s,site_tau_g,site_tau_f,site_tau_p)),
                    block = rnorm(4,0,c(block_tau_s,block_tau_g,block_tau_f,block_tau_p)),
                    source = rnorm(4,0,c(source_tau_s,source_tau_g,source_tau_f,source_tau_p)))
  rownames(rfx) <- c("surv","grow","flow","panic")
  return(rfx)
}


# set all parameters and run computations --------------------------------------

F_params <- M_params <- list()

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

# Estimation+process error
out <- lambdaSim_delay( F_params=F_params,
                        M_params=M_params,
                        grow_perturb=0,
                        surv_perturb=0,
                        flow_perturb=0,
                        fert_perturb=0,
                        viab_perturb=0,
                        zpptgrow=climate_values[l,5],
                        ztempgrow=climate_values[l,6],
                        zpptdorm=climate_values[l,7],
                        ztempdorm=climate_values[l,8],
                        rfx = rfx_fun(),
                        max.yrs = max_yrs)$lambdatracker[max_yrs]

# data frame to output
out_df <- climate_values[l,] %>% 
            mutate( lambda_diff  = out,
                    posterior_id = p )

# tell me how much time it took
print( Sys.time() - init_t )

# store result                      
write.csv( out_df, 
           paste0( '/gpfs1/work/compagna/poar_',
                   future_prj,'_',job_id,'.csv'), 
           row.names = F )
