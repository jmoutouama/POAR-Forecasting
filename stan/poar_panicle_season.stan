

data {
  // Data for all vital rates
  int<lower=0> n_sites;         // N. of sites
  int<lower=0> n_sources;         // N. of source populations

 

  // Data for fertility sub-model (p)
  int<lower=0> n_p;    // N. of data points for the fertility model
  int<lower=0> n_blocks_p;         // N. of blocks
  int<lower=0> site_p[n_p];  // site index
  int<lower=0> block_p[n_p];  // block index
  int<lower=0> source_p[n_p];  // source index
  int<lower=1> y_p[n_p]; //  number panicles at time t.
  vector[n_p] size_p;  // log size at time t
  vector[n_p] male_p;  // sex (male=1, fem=0)
  vector[n_p] pptgrow_p;  // precipitation of growing season 
  vector[n_p] pptdorm_p;  // precipitation of dormant season
  vector[n_p] tempgrow_p;  //  temperature of growing season
  vector[n_p] tempdorm_p;  //  temperature of growing season
}

parameters {
  //Panicles
  //fixed effects
  real b0_p;    
  real bsize_p;   
  real bsex_p;   
  real bpptgrow_p;  
  real bpptdorm_p;  
  real btempgrow_p;  
  real btempdorm_p;  
  real btempdormpptdorm_p;
  real btempgrowpptgrow_p;
  real bsizesex_p;
  real bpptgrowsex_p;  
  real bpptdormsex_p; 
  real btempgrowsex_p;  
  real btempdormsex_p;  
  real btempdormpptdormsex_p;  
  real btempgrowpptgrowsex_p;  

  //random effects
  real<lower=0> block_tau_p; 
  real block_rfx_p[n_blocks_p];  
  real<lower=0> source_tau_p; 
  real source_rfx_p[n_sources];
  real<lower=0> site_tau_p; 
  real site_rfx_p[n_sites];
  real<lower=0> phi_p; // Panicle dispersion parameter

  }

transformed parameters {

  real predP[n_p];
  
  
  // prediction for fertility
  for(ipan in 1:n_p){
    predP[ipan] = b0_p + 
                //main effects
                bsize_p * size_p[ipan] + bsex_p * male_p[ipan] + bpptgrow_p * pptgrow_p[ipan] + bpptdorm_p * pptdorm_p[ipan] + btempgrow_p * tempgrow_p[ipan]  + btempdorm_p * tempdorm_p [ipan] +
                
                //2-way interactions
                bpptgrowsex_p * pptgrow_p[ipan] * male_p[ipan] +
                bpptdormsex_p * pptdorm_p[ipan] * male_p[ipan] +
                btempgrowsex_p * tempgrow_p[ipan] * male_p[ipan] +
                btempdormsex_p * tempdorm_p[ipan] * male_p[ipan] +
                bsizesex_p * size_p[ipan] * male_p[ipan] +
                btempdormpptdorm_p * tempdorm_p[ipan] * pptdorm_p[ipan] +
                btempgrowpptgrow_p * tempgrow_p[ipan] * tempgrow_p[ipan] +

                //3-way interaction
                btempdormpptdormsex_p * tempdorm_p[ipan] * pptdorm_p[ipan] * male_p[ipan] +
                btempgrowpptgrowsex_p * tempgrow_p[ipan] * tempgrow_p[ipan] * male_p[ipan] +


                //random effects
                block_rfx_p[block_p[ipan]] +
                source_rfx_p[source_p[ipan]]+
                site_rfx_p[site_p[ipan]];
    }

}

model {
  // priors on parameters
  // fertility
  b0_p ~ normal(0, 500);    
  bsize_p ~ normal(0, 100);   
  bsizesex_p ~ normal(0, 100);   
  bsex_p ~ normal(0, 100);   
  bpptgrow_p ~ normal(0, 100);  
  bpptdorm_p ~ normal(0, 100);  
  btempgrow_p ~ normal(0, 100);  
  btempdorm_p ~ normal(0, 100);
  bpptgrowsex_p ~ normal(0, 100);
  bpptdormsex_p ~ normal(0, 100);
  btempgrowsex_p ~ normal(0, 100);
  btempdormsex_p ~ normal(0, 100);
  btempdormpptdorm_p ~ normal(0, 100);
  btempgrowpptgrow_p ~ normal(0, 100);
  btempdormpptdormsex_p ~ normal(0, 100);
  btempgrowpptgrowsex_p ~ normal(0, 100);
  block_tau_p ~ inv_gamma(0.2, 0.2);
  for (i in 1:n_blocks_p){
    block_rfx_p[i] ~ normal(0, block_tau_p);
  }
  source_tau_p ~ inv_gamma(0.2, 0.2);
  for (i in 1:n_sources){
    source_rfx_p[i] ~ normal(0, source_tau_p);
  }
  site_tau_p ~ inv_gamma(0.2, 0.2);
  for (i in 1:n_sites){
    site_rfx_p[i] ~ normal(0, site_tau_p);
  }

  // sampling  

  //fertility need loop for zero truncation
  for (i in 1:n_p) {
    y_p[i] ~ neg_binomial_2_log(predP[i], phi_p);
    // manually zero-truncating
    target += - log1m(neg_binomial_2_log_lpmf(0 | predP[i], phi_p)); 
  }
  
}

generated quantities {
  vector[n_p] log_lik;
  for (i in 1:n_p) {
  log_lik[i] = neg_binomial_2_log_lpmf (y_p[i] | predP[i],phi_p);
 }
}



