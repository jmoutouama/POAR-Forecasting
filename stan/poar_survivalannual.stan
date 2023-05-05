
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
  vector[n_s] ppt_s;  // precipitation of site
  vector[n_s] temp_s;  //  temperature of site
  vector[n_s] cvppt_s;  // precipitation coefficient of variation of site
  
}

parameters {
   //Survival
  //fixed effects
  real b0_s;    
  real bsize_s;   
  real bsex_s;   
  real bppt_s;  
  real btemp_s;  
  real bcvppt_s;  
  real bpptsex_s;
  real btempsex_s;
  real bcvpptsex_s;
  real bppttempsex_s;
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
                bsize_s * size_s[isurv] + bsex_s * male_s[isurv] + bppt_s * ppt_s[isurv] + btemp_s * temp_s[isurv]+ bcvppt_s * cvppt_s[isurv] +
                //2-way interactions
                bpptsex_s * ppt_s[isurv] * male_s[isurv] +
                btempsex_s * temp_s[isurv] * male_s[isurv] +
                bcvpptsex_s * cvppt_s[isurv] * male_s[isurv] +
                
                //3-way interaction
                bppttempsex_s * temp_s[isurv] * ppt_s[isurv] * male_s[isurv] +
               

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
  bppt_s ~ normal(0, 100);  
  btemp_s ~ normal(0, 100);  
  bcvppt_s ~ normal(0, 100);  
  bpptsex_s ~ normal(0, 100);
  btempsex_s ~ normal(0, 100);
  bppttempsex_s ~ normal(0, 100);
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
 log_lik[i] = bernoulli_logit_lpmf(y_s[i] | predS[i]);
 }
}

   



