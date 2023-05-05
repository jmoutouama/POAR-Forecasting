
data {
  // Data for all vital rates
  int<lower=0> n_sites;         // N. of sites
  int<lower=0> n_sources;         // N. of source populations

  // Data for flowering sub-model (f)
  int<lower=0> n_f;    // N. of data points for the flowering  model
  int<lower=0> n_blocks_f;         // N. of blocks
  int<lower=0> site_f[n_f];  // site index
  int<lower=0> block_f[n_f];  // block index
  int<lower=0> source_f[n_f];  // source index
  int<lower=0,upper=1> y_f[n_f]; // Flowering at time t+1.
  vector[n_f] size_f;  // log size at time t
  vector[n_f] male_f;  // sex (male=1, fem=0)
  vector[n_f] ppt_f;  // precipitation of site
  vector[n_f] temp_f;  //  temperature of site
  vector[n_f] cvppt_f;  // precipitation standard deviation  of site

}

parameters {
  //Flowering
  //fixed effects
  real b0_f;    
  real bsize_f;   
  real bsex_f;   
  real bppt_f;  
  real btemp_f;  
  real bcvppt_f;  
  real bpptsex_f;
  real btempsex_f;
  real bcvpptsex_f;
  real bppttempsex_f;
  //random effects
  real<lower=0> block_tau_f; 
  real block_rfx_f[n_blocks_f];  
  real<lower=0> source_tau_f; 
  real source_rfx_f[n_sources];
  real<lower=0> site_tau_f; 
  real site_rfx_f[n_sites];

  }

transformed parameters {
  real predF[n_f];
  // prediction for flowering
  for(iflow in 1:n_f){
    predF[iflow] = b0_f + 
                //main effects
                bsize_f * size_f[iflow] + bsex_f * male_f[iflow] + bppt_f * ppt_f[iflow] + btemp_f * temp_f[iflow]+ bcvppt_f * cvppt_f[iflow] + 
                //2-way interactions
                bpptsex_f * ppt_f[iflow] * male_f[iflow] +
                btempsex_f * temp_f[iflow] * male_f[iflow] +
                bcvpptsex_f * cvppt_f[iflow] * male_f[iflow] +
               
                //3-way interaction
                bppttempsex_f * temp_f[iflow] * ppt_f[iflow] * male_f[iflow] +
             

                //random effects
                block_rfx_f[block_f[iflow]] +
                source_rfx_f[source_f[iflow]]+
                site_rfx_f[site_f[iflow]];
    }
}

model {
  // priors on parameters
  
  
  // Flowering
  b0_f ~ normal(0, 500);    
  bsize_f ~ normal(0, 100);   
  bsex_f ~ normal(0, 100);   
  bppt_f ~ normal(0, 100);  
  btemp_f ~ normal(0, 100);  
  bcvppt_f ~ normal(0, 100);  
  bpptsex_f ~ normal(0, 100);
  btempsex_f ~ normal(0, 100);
  bppttempsex_f ~ normal(0, 100);
  block_tau_f ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_blocks_f){
    block_rfx_f[i] ~ normal(0, block_tau_f);
  }
  source_tau_f ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_sources){
    source_rfx_f[i] ~ normal(0, source_tau_f);
  }
  site_tau_f ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_sites){
    site_rfx_f[i] ~ normal(0, site_tau_f);
  }

  // sampling  

  
  //flowering
  y_f ~ bernoulli_logit(predF);
  
}

generated quantities {
  vector[n_f] log_lik;
  for (i in 1:n_f) {
 log_lik[i] = bernoulli_logit_lpmf(y_f[i] | predF[i]);
 }
}


