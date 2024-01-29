
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
  grow<-dpoisinvgauss(x=y,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  grow<-ifelse(is.nan(grow) | is.infinite(grow),0,grow)
  truncLower<-dpoisinvgauss(x=0,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  truncUpper<-dpoisinvgauss(x=params$max_size:10000,mean=grow_mean,shape=(grow_mean*params$sigma_g))
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

