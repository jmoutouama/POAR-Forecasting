
# vital rate and megamatrix functions ---------------------------------------------------
invlogit<-function(x){exp(x)/(1+exp(x))}

#SURVIVAL AT SIZE X.
sx<-function(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,surv_perturb=0){
  # x_scaled_s=(x-7.768868)/9.855441
  pptgrow_scaled_s=(pptgrow-752.6118)/336.1784
  pptdorm_scaled_s=(pptdorm-379.6936)/146.842
  tempgrow_scaled_s=(tempgrow- 15.02675)/3.137914
  tempdorm_scaled_s=(tempdorm-27.36059)/3.137914
  surv_mean<-params$surv_mu + 
    params$surv_size*log(x) + 
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
  return(invlogit(surv_mean)+ surv_perturb)
}

# plot(sx(x=1:40,
#    params=M_params,
#    pptgrow,pptdorm,tempgrow,tempdorm,rfx))


#PROBABILITY OF GROWTH FROM SIZE X TO Y
#This function truncates the density asscociation with x==0 and x>x.max
gxy<-function(x,y,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,grow_perturb=0){
    #x_scaled_g=(x-9.019255)/10.83937
    pptgrow_scaled_g=(pptgrow-709.7987)/338.7739
    pptdorm_scaled_g=(pptdorm-369.8298)/143.5237
    tempgrow_scaled_g=(tempgrow-15.16865)/3.149584
    tempdorm_scaled_g=(tempdorm- 27.34114)/1.468994
  grow_mean<-exp(params$grow_mu + 
    params$grow_size*log(x) + 
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
    rfx["grow","site"] + rfx["grow","block"] + rfx["grow","source"]) + grow_perturb
  # grow<-dpoisinvgauss(x=y,mean=grow_mean,shape=grow_mean*params$sigma_g)
  # truncZero<-dpoisinvgauss(x=0,mean=grow_mean,shape=grow_mean*params$sigma_g)
  # return(grow/(1-truncZero))
  grow<-dpoisinvgauss(x=y,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  grow<-ifelse(is.nan(grow) | is.infinite(grow),0,grow)
  truncLower<-dpoisinvgauss(x=0,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  truncUpper<-dpoisinvgauss(x=params$max_size:10000,mean=grow_mean,shape=(grow_mean*params$sigma_g))
  truncUpper<-sum(ifelse(is.nan(truncUpper) | is.infinite(truncUpper),0,truncUpper))
  return(grow/(1-(truncLower+truncUpper)))
}


 
# sum(gxy(100,1:120,params=F_params,pptgrow,pptdorm,tempgrow,tempdorm,rfx))
# sum(gxy(F_params$max_size,1:1000,params=F_params,pptgrow,pptdorm,tempgrow,tempdorm,rfx))
# y=1:40
# plot(gxy(x=y,y=y,params=F_params,pptgrow,pptdorm,tempgrow,tempdorm,rfx=rfx,grow_perturb=0))

#SURVIVAL*GROWTH 
pxy<-function(x,y,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,surv_perturb=0,grow_perturb=0){
  sx(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,surv_perturb) * gxy(x,y,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,grow_perturb)
}

# y=1:40
# plot(pxy(x=y,y=y,params=F_params,pptgrow,pptdorm,tempgrow,tempdorm,rfx=rfx,surv_perturb=0))

# PROBABILITY OF FLOWERING
pfx<-function(x,params, pptgrow,pptdorm,tempgrow,tempdorm,rfx,flow_perturb=0){
  # x_scaled_f=(x-19.39529)/38.72338
  pptgrow_scaled_f=(pptgrow-706.8709 )/339.4783
  pptdorm_scaled_f=(pptdorm-366.5321)/140.0161
  tempgrow_scaled_f=(tempgrow-15.2168)/3.142949
  tempdorm_scaled_f=(tempdorm-27.34905)/1.453832
  flow_mean<-params$flow_mu + 
    params$flow_size*log(x) + 
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
  return(invlogit(flow_mean) + flow_perturb )
}

# y=1:40
# plot(pfx(x=y,params=F_params,pptgrow,pptdorm,tempgrow,tempdorm,rfx=rfx,flow_perturb=0))

#NUMBER OF PANICLES
nfx<-function(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,fert_perturb=0){
  # x_scaled_p=(x-38.78495)/58.04579
  pptgrow_scaled_p=(pptgrow-564.0412)/273.6934 
  pptdorm_scaled_p=(pptdorm-342.8549)/134.8506
  tempgrow_scaled_p=(tempgrow-14.40408)/2.964567
  tempdorm_scaled_p=(tempdorm-27.02796)/1.43082
  panic_mean<-params$panic_mu + 
    params$panic_size*log(x) + 
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
  return(exp(panic_mean) + fert_perturb)
}

# y=1:40
# plot(pfx(x=y,params=F_params,pptgrow,pptdorm,tempgrow,tempdorm,rfx=rfx,flow_perturb=0))

#SEED VIABILITY
viab<-function(params,twosex,OSR=NULL,viab_perturb=0){
  if(twosex==F){return(params$v0 + viab_perturb)}
  if(twosex==T){return((params$v0 * (1 - OSR ^ params$a_v))+ viab_perturb)}
}

#FERTILITY--returns number of recruits
##Female offspring
fertx_F<-function(x,params,rfx,pptgrow,pptdorm,tempgrow,tempdorm,twosex,OSR=NULL,flow_perturb=0,fert_perturb=0,viab_perturb=0){
  seedlings<-pfx(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,flow_perturb)*nfx(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,fert_perturb)*params$ov_per_inf*viab(params,twosex,OSR,viab_perturb)*params$germ*params$PSR
  return(seedlings)
}

# y=1:40
# plot(fertx_F(x=y,params=F_params,rfx=rfx,twosex=F,OSR=NULL,pptgrow,pptdorm,tempgrow,tempdorm,flow_perturb=0,fert_perturb=0,viab_perturb=0))


##Male offspring
fertx_M<-function(x,params,rfx,pptgrow,pptdorm,tempgrow,tempdorm,twosex,OSR=NULL,flow_perturb=0,fert_perturb=0,viab_perturb=0){
  seedlings<-pfx(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,flow_perturb)*nfx(x,params,pptgrow,pptdorm,tempgrow,tempdorm,rfx,fert_perturb)*params$ov_per_inf*viab(params,twosex,OSR,viab_perturb)*params$germ*(1-params$PSR)
  return(seedlings)
}

# y=1:40
# plot(fertx_M(x=y,params=F_params,rfx=rfx,twosex=F,OSR=NULL,pptgrow,pptdorm,tempgrow,tempdorm,flow_perturb=0,fert_perturb=0,viab_perturb=0))

megamatrix_delay<-function(F_params,M_params,pptgrow,pptdorm,tempgrow,tempdorm,twosex,OSR=NULL,rfx,
                           surv_perturb=0,grow_perturb=0,flow_perturb=0,fert_perturb=0,viab_perturb=0){  
  matdim<-F_params$max_size       
  y<-1:F_params$max_size
  
  ## F-to-F (growth/survival transition)
  F.Tmat<-matrix(0,matdim+1,matdim+1)
  F.Tmat[2:(matdim+1),2:(matdim+1)]<-t(outer(y,y,pxy,params=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx,surv_perturb=surv_perturb, grow_perturb=grow_perturb))
  # surviving seedlings emerge in continuous population
  F.Tmat[2:(matdim+1),1] <- gxy(x=1,y=y,params=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx,grow_perturb=grow_perturb) * (M_params$sdlg_surv + surv_perturb)
  
  # F-to-F Fertility transition
  F.Fmat<-matrix(0,matdim+1,matdim+1)
  # seedlings in top row
  F.Fmat[1,2:(matdim+1)]<-fertx_F(x=y,params=F_params,rfx=rfx,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,twosex=twosex,OSR=OSR,flow_perturb=flow_perturb,fert_perturb=fert_perturb,viab_perturb=viab_perturb) 
  
  ## M-to-M (growth/survival transition)
  M.Tmat<-matrix(0,matdim+1,matdim+1)
  M.Tmat[2:(matdim+1),2:(matdim+1)]<-t(outer(y,y,pxy,params=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx,surv_perturb=surv_perturb, grow_perturb=grow_perturb))
  M.Tmat[2:(matdim+1),1] <- gxy(x=1,y=y,params=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx,grow_perturb=grow_perturb) * (M_params$sdlg_surv + surv_perturb)
  
  # F-to-M Fertility transition
  M.Fmat<-matrix(0,matdim+1,matdim+1)
  M.Fmat[1,2:(matdim+1)]<-fertx_M(x=y,params=F_params,rfx=rfx,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,twosex=twosex,OSR=OSR,flow_perturb=flow_perturb,fert_perturb=fert_perturb,viab_perturb=viab_perturb) 
  
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
                          grow_perturb=0,surv_perturb=0,flow_perturb=0,fert_perturb=0,viab_perturb=0){
  matdim<-F_params$max_size       
  y<-1:F_params$max_size
  lambdatracker      <- rep(0,max.yrs)
  OSRtracker   <- rep(0,max.yrs)
  SRtracker   <- rep(0,max.yrs)
  n0            <- rep(1/((matdim+1)*2),((matdim+1)*2))
  
  for(t in 1:max.yrs){
    ##Estimate panicle SR
    flowering_females<-n0[2:(matdim+1)]*pfx(x=y,param=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx,flow_perturb=flow_perturb) ## scalar multiplication to weight females by flowering prob
    F_panicles<-flowering_females%*%nfx(x=y,param=F_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx,fert_perturb=fert_perturb) ##Vector operation to sum female panicles
    flowering_males<-n0[(matdim+3):((matdim+1)*2)]*pfx(x=y,param=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx,flow_perturb=flow_perturb)
    M_panicles<-flowering_males %*% nfx(x=y,param=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,rfx=rfx,fert_perturb=fert_perturb)
    OSRtracker[t]<-F_panicles/(F_panicles+M_panicles) ##Panicle sex ratio (proportion female)
    SRtracker[t]<-sum(n0[1:(matdim+1)])
    #assmble matrix
    MEGAmat<-megamatrix_delay(F_params=F_params,M_params=M_params,pptgrow=pptgrow,pptdorm=pptdorm,tempgrow=tempgrow,tempdorm=tempdorm,twosex=T,OSR=OSRtracker[t],rfx=rfx,
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

