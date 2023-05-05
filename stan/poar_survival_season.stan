data {
  // Data for all vital rates
  int<lower=0> n_sites;         // N. of sites
  int<lower=0> n_sources;         // N. of source populations

  // Data for survival sub-model (s)
  int<lower=0> n_s;    // N. of data points for the surival  model
  int<lower=0> n_blocks_s;         // N. of blocks
  int<lower=0> site_s[n_s];  // site index
  int<lower=0> block_s[n_s];  // block index
  int<lower=0> source_s[n_s];  // source index
  int<lower=0,upper=1> y_s[n_s]; // Survival at time t+1.
  vector[n_s] size_s;  // log size at time t
  vector[n_s] male_s;  // sex (male=1, fem=0)
  vector[n_s] pptgrow_s;  // precipitation of growing season 
  vector[n_s] pptdorm_s;  // precipitation of dormant season
  vector[n_s] tempgrow_s;  //  temperature of growing season
  vector[n_s] tempdorm_s;  //  temperature of growing season
  
}

parameters {
   //Survival
  //fixed effects
  real b0_s;    
  real bsize_s;   
  real bsex_s;   
  real bpptgrow_s;  
  real bpptdorm_s;  
  real btempgrow_s;  
  real btempdorm_s;  
  real bpptgrowsex_s;  
  real bpptdormsex_s; 
  real btempgrowsex_s;  
  real btempdormsex_s;  
  real btempdormpptdormsex_s;  
  real btempgrowpptgrowsex_s;  

  //random effects
  real<lower=0> block_tau_s; 
  real block_rfx_s[n_blocks_s];  
  real<lower=0> source_tau_s; 
  real source_rfx_s[n_sources];
  real<lower=0> site_tau_s; 
  real site_rfx_s[n_sites];
  }

transformed parameters {

  real predS[n_s];

  // prediction for survival
  for(isurv in 1:n_s){
    predS[isurv] = b0_s + 
                //main effects
                bsize_s * size_s[isurv] + bsex_s * male_s[isurv] + bpptgrow_s * pptgrow_s[isurv] + bpptdorm_s * pptdorm_s[isurv] + btempgrow_s * tempgrow_s[isurv]  + btempdorm_s * tempdorm_s [isurv] +
                
                //2-way interactions
                bpptgrowsex_s * pptgrow_s[isurv] * male_s[isurv] +
                bpptdormsex_s * pptdorm_s[isurv] * male_s[isurv] +
                btempgrowsex_s * tempgrow_s[isurv] * male_s[isurv] +
                btempdormsex_s * tempdorm_s[isurv] * male_s[isurv] +

                //3-way interaction
                btempdormpptdormsex_s * tempdorm_s[isurv] * pptdorm_s[isurv] * male_s[isurv] +
                btempgrowpptgrowsex_s * tempgrow_s[isurv] * tempgrow_s[isurv] * male_s[isurv] +
                
                //random effects
                block_rfx_s[block_s[isurv]] +
                source_rfx_s[source_s[isurv]]+
                site_rfx_s[site_s[isurv]];
    }

}

model {
  // priors on parameters
  //Survival
  b0_s ~ normal(0, 100);    
  bsize_s ~ normal(0, 100);   
  bsex_s ~ normal(0, 100);   
  bpptgrow_s ~ normal(0, 100);  
  bpptdorm_s ~ normal(0, 100);  
  btempgrow_s ~ normal(0, 100);  
  btempdorm_s ~ normal(0, 100);
  bpptgrowsex_s ~ normal(0, 100);
  bpptdormsex_s ~ normal(0, 100);
  btempgrowsex_s ~ normal(0, 100);
  btempdormsex_s ~ normal(0, 100);
  block_tau_s ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_blocks_s){
    block_rfx_s[i] ~ normal(0, block_tau_s);
  }
  source_tau_s ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_sources){
    source_rfx_s[i] ~ normal(0, source_tau_s);
  }
  site_tau_s ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_sites){
    site_rfx_s[i] ~ normal(0, site_tau_s);
  }

  // sampling  

  //survival
  y_s ~ bernoulli_logit(predS);
}

generated quantities {
  vector[n_s] log_lik;
  for (i in 1:n_s) {
 log_lik[i] = bernoulli_logit_lpmf(y_s[i] |predS[i]);
 }
}




