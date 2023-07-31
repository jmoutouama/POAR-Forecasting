
# vital rate and megamatrix functions ---------------------------------------------------

#SURVIVAL AT SIZE X.
sx<-function(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx){
  x_scaled_s=(x-7.768868)/9.855441
  pptgrow_scaled_s=(pptgrow-1425.061)/536.1011
  pptdorm_scaled_s=(pptdorm-755.4659)/240.2107
  tempgrow_scaled_s=(tempgrow-15.12044)/3.139186
  tempdorm_scaled_s=(tempdorm-27.35746)/1.396525
  surv_mean<-params$surv_mu + 
    params$surv_size*x_scaled_s + 
    params$surv_pptgrow*pptgrow_scaled_s + 
    params$surv_pptdorm*pptdorm_scaled_s + 
    params$surv_tempgrow*tempgrow_scaled_s + 
    params$surv_tempdorm*tempdorm_scaled_s + 
    params$surv_tempdorm_pptdorm*tempdorm_scaled_s*pptdorm_scaled_s + 
    params$surv_tempgrow_pptgrow*tempgrow_scaled_s*pptgrow_scaled_s + 
    params$surv_pptgrow2*(pptgrow_scaled_s^2) + 
    params$surv_pptdorm2*(pptdorm_scaled_s^2) + 
    params$surv_tempgrow2*(tempgrow_scaled_s^2) + 
    params$surv_tempdorm2*(tempdorm_scaled_s^2) +
    rfx["surv","site"] + rfx["surv","block"] + rfx["surv","source"]
  return(invlogit(surv_mean))
}

# y=1:40
# plot(sx(x=y,params=F_params,pptgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptgrow)),pptdorm=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptdorm)),tempgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$tempgrow)),tempdorm=mean(poarclim_seasonal.sur_mean_sd$tempdorm),rfx=rfx,surv_perturb=0))


#PROBABILITY OF GROWTH FROM SIZE X TO Y
#This function truncates the density asscociation with x==0 and x>x.max
gxy<-function(x,y,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx){
    x_scaled_g=(x-9.019255)/10.83937
    pptgrow_scaled_g=(pptgrow-1384.342)/562.9291
    pptdorm_scaled_g=(pptdorm-734.1461)/232.1593
    tempgrow_scaled_g=(tempgrow-15.21964)/3.133994
    tempdorm_scaled_g=(tempdorm-27.33936)/1.433536
  grow_mean<-exp(params$grow_mu + 
    params$grow_size*x_scaled_g + 
    params$grow_pptgrow*pptgrow_scaled_g + 
    params$grow_pptdorm*pptdorm_scaled_g + 
    params$grow_tempgrow*tempgrow_scaled_g + 
    params$grow_tempdorm*tempdorm_scaled_g + 
    params$grow_tempdorm_pptdorm*tempdorm_scaled_g*pptdorm_scaled_g + 
    params$grow_tempgrow_pptgrow*tempgrow_scaled_g*pptgrow_scaled_g + 
    params$grow_pptgrow2*(pptgrow_scaled_g^2) + 
    params$grow_pptdorm2*(pptdorm_scaled_g^2) + 
    params$grow_tempgrow2*(tempgrow_scaled_g^2) + 
    params$grow_tempdorm2*(tempdorm_scaled_g^2) +
    rfx["grow","site"] + rfx["grow","block"] + rfx["grow","source"]) 
  
  grow<-dpoisinvgauss(x=y,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  grow<-ifelse(is.nan(grow) | is.infinite(grow),0,grow)
  truncLower<-dpoisinvgauss(x=0,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  truncUpper<-dpoisinvgauss(x=params$max_size:10000,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  truncUpper<-sum(ifelse(is.nan(truncUpper) | is.infinite(truncUpper),0,truncUpper))
  return(grow/(1-(truncLower+truncUpper)))
}


# y=1:40
# plot(gxy(x=y,y=y,params=F_params,pptgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptgrow)),pptdorm=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptdorm)),tempgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$tempgrow)),tempdorm=mean(poarclim_seasonal.sur_mean_sd$tempdorm),rfx=rfx,grow_perturb=0))

#SURVIVAL*GROWTH 
pxy<-function(x,y,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx){
  sx(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx) * gxy(x,y,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx)
}

# y=1:40
# plot(pxy(x=y,y=y,params=F_params,pptgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptgrow)),pptdorm=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptdorm)),tempgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$tempgrow)),tempdorm=mean(poarclim_seasonal.sur_mean_sd$tempdorm),rfx=rfx,surv_perturb=0))



# PROBABILITY OF FLOWERING
pfx<-function(x,params, pptgrow,pptdorm,tempgrow,tempdorm,rfx){
  x_scaled_f=(x-19.39529)/38.72338
  pptgrow_scaled_f=(pptgrow-1388.22)/567.5876
  pptdorm_scaled_f=(pptdorm-731.769)/228.7285
  tempgrow_scaled_f=(tempgrow-15.25128)/3.12259
  tempdorm_scaled_f=(tempdorm-27.34895)/1.426207
  flow_mean<-params$flow_mu + 
    params$flow_size*x_scaled_f + 
    params$flow_pptgrow*pptgrow_scaled_f + 
    params$flow_pptdorm*pptdorm_scaled_f + 
    params$flow_tempgrow*tempgrow_scaled_f + 
    params$flow_tempdorm*tempdorm_scaled_f + 
    params$flow_tempdorm_pptdorm*tempdorm_scaled_f*pptdorm_scaled_f + 
    params$flow_tempgrow_pptgrow*tempgrow_scaled_f*pptgrow_scaled_f + 
    params$flow_pptgrow2*(pptgrow_scaled_f^2) + 
    params$flow_pptdorm2*(pptdorm_scaled_f^2) + 
    params$flow_tempgrow2*(tempgrow_scaled_f^2) + 
    params$flow_tempdorm2*(tempdorm_scaled_f^2) +
    rfx["flow","site"] + rfx["flow","block"] + rfx["flow","source"]
  return(invlogit(flow_mean) )
}

# y=1:40
# plot(pfx(x=y,params=F_params,pptgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptgrow)),pptdorm=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptdorm)),tempgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$tempgrow)),tempdorm=mean(poarclim_seasonal.sur_mean_sd$tempdorm),rfx=rfx,flow_perturb=0))

#NUMBER OF PANICLES
nfx<-function(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx){
  x_scaled_p=(x-38.78495)/58.04579
  pptgrow_scaled_p=(pptgrow-1119.937)/440.4973
  pptdorm_scaled_p=(pptdorm-679.8125)/221.7003
  tempgrow_scaled_p=(tempgrow-14.41257)/2.941214
  tempdorm_scaled_p=(tempdorm-27.04046)/1.42038
  panic_mean<-params$panic_mu + 
    params$panic_size*x_scaled_p + 
    params$panic_pptgrow*pptgrow_scaled_p + 
    params$panic_pptdorm*pptdorm_scaled_p + 
    params$panic_tempgrow*tempgrow_scaled_p + 
    params$panic_tempdorm*tempdorm_scaled_p + 
    params$panic_tempdorm_pptdorm*tempdorm_scaled_p*pptdorm_scaled_p + 
    params$panic_tempgrow_pptgrow*tempgrow_scaled_p*pptgrow_scaled_p + 
    params$panic_pptgrow2*(pptgrow_scaled_p^2) + 
    params$panic_pptdorm2*(pptdorm_scaled_p^2) + 
    params$panic_tempgrow2*(tempgrow_scaled_p^2) + 
    params$panic_tempdorm2*(tempdorm_scaled_p^2) +
    rfx["panic","site"] + rfx["panic","block"] + rfx["panic","source"]
  return(exp(panic_mean) )
}

# y=1:40
# plot(pfx(x=y,params=F_params,pptgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptgrow)),pptdorm=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptdorm)),tempgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$tempgrow)),tempdorm=mean(poarclim_seasonal.sur_mean_sd$tempdorm),rfx=rfx,flow_perturb=0))

#SEED VIABILITY
viab<-function(params,twosex,OSR=NULL){
  if(twosex==F){return(params$v0 )}
  if(twosex==T){return((params$v0 * (1 - OSR ^ params$a_v)))}
}

#FERTILITY--returns number of recruits
##Female offspring
fertx_F<-function(x,params,rfx,pptgrow,pptdorm,tempgrow,tempdorm,twosex,OSR=NULL){
  seedlings<-pfx(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx)*nfx(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx)*params$ov_per_inf*viab(params,twosex,OSR)*params$germ*params$PSR
  return(seedlings)
}

# y=1:40
# plot(fertx_F(x=y,params=F_params,rfx=rfx,twosex=F,OSR=NULL,pptgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptgrow)),pptdorm=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptdorm)),tempgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$tempgrow)),tempdorm=mean(poarclim_seasonal.sur_mean_sd$tempdorm)))
# # 

##Male offspring
fertx_M<-function(x,params,rfx,pptgrow,pptdorm,tempgrow,tempdorm,twosex,OSR=NULL){
  seedlings<-pfx(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx)*nfx(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx)*params$ov_per_inf*viab(params,twosex,OSR)*params$germ*(1-params$PSR)
  return(seedlings)
}

# y=1:40
# plot(fertx_M(x=y,params=M_params,rfx=rfx,twosex=F,OSR=NULL,pptgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptgrow)),pptdorm=mean(na.omit(poarclim_seasonal.sur_mean_sd$pptdorm)),tempgrow=mean(na.omit(poarclim_seasonal.sur_mean_sd$tempgrow)),tempdorm=mean(poarclim_seasonal.sur_mean_sd$tempdorm),flow_perturb=0,fert_perturb=0,viab_perturb=0))

megamatrix_delay<-function(F_params,M_params,pptgrow,pptdorm,tempgrow,tempdorm,twosex,OSR=NULL,rfx){  
  matdim<-F_params$max_size       
  y<-1:F_params$max_size
  
  ## F-to-F (growth/survival transition)
  F.Tmat<-matrix(0,matdim+1,matdim+1)
  F.Tmat[2:(matdim+1),2:(matdim+1)]<-t(outer(y,y,pxy,params=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx))
  # surviving seedlings emerge in continuous population
  F.Tmat[2:(matdim+1),1] <- gxy(x=1,y=y,params=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx) * (M_params$sdlg_surv)
  
  # F-to-F Fertility transition
  F.Fmat<-matrix(0,matdim+1,matdim+1)
  # seedlings in top row
  F.Fmat[1,2:(matdim+1)]<-fertx_F(x=y,params=F_params,rfx=rfx,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,twosex=twosex,OSR=OSR) 
  
  ## M-to-M (growth/survival transition)
  M.Tmat<-matrix(0,matdim+1,matdim+1)
  M.Tmat[2:(matdim+1),2:(matdim+1)]<-t(outer(y,y,pxy,params=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx))
  M.Tmat[2:(matdim+1),1] <- gxy(x=1,y=y,params=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx) * (M_params$sdlg_surv)
  
  # F-to-M Fertility transition
  M.Fmat<-matrix(0,matdim+1,matdim+1)
  M.Fmat[1,2:(matdim+1)]<-fertx_M(x=y,params=F_params,rfx=rfx,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,twosex=twosex,OSR=OSR) 
  
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

lambdaSim_delay<-function(F_params,M_params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,max.yrs,
                          grow_perturb=0,surv_perturb=0,flow_perturb=0,fert_perturb=0){
  matdim<-F_params$max_size       
  y<-1:F_params$max_size
  lambdatracker      <- rep(0,max.yrs)
  OSRtracker   <- rep(0,max.yrs)
  SRtracker   <- rep(0,max.yrs)
  n0            <- rep(1/((matdim+1)*2),((matdim+1)*2))
  
  for(t in 1:max.yrs){
    ##Estimate panicle SR
    flowering_females<-n0[2:(matdim+1)]*pfx(x=y,param=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx) ## scalar multiplication to weight females by flowering prob
    F_panicles<-flowering_females%*%nfx(x=y,param=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx) ##Vector operation to sum female panicles
    flowering_males<-n0[(matdim+3):((matdim+1)*2)]*pfx(x=y,param=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx)
    M_panicles<-flowering_males%*%nfx(x=y,param=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx)
    OSRtracker[t]<-F_panicles/(F_panicles+M_panicles) ##Panicle sex ratio (proportion female)
    SRtracker[t]<-sum(n0[1:(matdim+1)])
    #assmble matrix
    MEGAmat<-megamatrix_delay(F_params=F_params,M_params=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,twosex=T,OSR=OSRtracker[t],rfx=rfx
                              )$MEGAmat
    n0 <- MEGAmat[,] %*% n0
    N  <- sum(n0)
    lambdatracker[t]<-N
    n0 <-n0/N
  }
  return(list(lambdatracker=lambdatracker,SRtracker=SRtracker,OSRtracker=OSRtracker,n0=n0))
}



# RFX fun -----------------------------------------------------------------
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

# basement garbage --------------------------------------------------------


## assmble vital rate functions into a projection matrix
## This was the original megamatrix function, was not updated to include vital rate perturbations (because we are not using it)
megamatrix<-function(F_params,M_params,pptgrow,pptdorm,tempgrow,tempdorm,twosex,OSR=NULL,rfx){  
  matdim<-F_params$max_size         
  y<-1:F_params$max_size
  
  ## F-to-F (growth/survival transition)
  F.Tmat<-matrix(0,matdim,matdim)
  F.Tmat[1:matdim,1:matdim]<-t(outer(y,y,pxy,params=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx))
  
  # F-to-F Fertility transition
  F.Fmat<-matrix(0,matdim,matdim)
  F.Fmat[1,1:matdim]<-fertx_F(x=y,params=F_params,rfx=rfx,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,twosex=twosex,OSR=OSR) 
  
  ## M-to-M (growth/survival transition)
  M.Tmat<-matrix(0,matdim,matdim)
  M.Tmat[1:matdim,1:matdim]<-t(outer(y,y,pxy,params=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx))
  
  # F-to-M Fertility transition
  M.Fmat<-matrix(0,matdim,matdim)
  M.Fmat[1,1:matdim]<-fertx_M(x=y,params=F_params,rfx=rfx,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,twosex=twosex,OSR=OSR) 
  
  #M-to-F
  zero.mat<-matrix(0,matdim,matdim)
  
  # Put it all together as a megamatrix
  MEGAmat<-cbind(rbind(F.Tmat+F.Fmat,  ##Female growth/survival + recruitment[1,1]
                       M.Fmat), ##Male recruitment [2,1]
                 rbind(zero.mat,   ##Females from males [1,2]
                       M.Tmat))   ##Male growth/survival
  
  return(list(MEGAmat=MEGAmat,F.Tmat=F.Tmat,M.Tmat=M.Tmat,y=y))
}


lambdaSim<-function(F_params,M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx,max.yrs){
  matdim<-F_params$max_size         
  y<-1:F_params$max_size
  lambdatracker      <- rep(0,max.yrs)
  OSRtracker   <- rep(0,max.yrs)
  SRtracker   <- rep(0,max.yrs)
  n0            <- rep(1/(matdim*2),(matdim*2))
  
  for(t in 1:max.yrs){
    ##Estimate panicle SR
    flowering_females<-n0[1:matdim]*pfx(x=y,param=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx) ## scalar multiplication to weight females by flowering prob
    F_panicles<-flowering_females%*%nfx(x=y,param=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx) ##Vector operation to sum female panicles
    flowering_males<-n0[(matdim+1):(matdim*2)]*pfx(x=y,param=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx)
    M_panicles<-flowering_males%*%nfx(x=y,param=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx)
    OSRtracker[t]<-F_panicles/(F_panicles+M_panicles) ##Panicle sex ratio (proportion female)
    SRtracker[t]<-sum(n0[1:matdim])
    #assmble matrix
    MEGAmat<-megamatrix(F_params=F_params,M_params=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,twosex=T,OSR=OSRtracker[t],rfx=rfx)$MEGAmat
    n0 <- MEGAmat[,] %*% n0
    N  <- sum(n0)
    lambdatracker[t]<-N
    n0 <-n0/N
  }
  return(list(lambdatracker=lambdatracker,SRtracker=SRtracker,OSRtracker=OSRtracker,n0=n0))
}


