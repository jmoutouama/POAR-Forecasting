

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
  int<lower=1> y_p[n_p]; // # panicles at time t.
  vector[n_p] size_p;  // log size at time t
  vector[n_p] male_p;  // sex (male=1, fem=0)
  vector[n_p] ppt_p;  // precipitation of site
  vector[n_p] temp_p;  //  temperature of site
  vector[n_p] cvppt_p;  // precipitation coefficient of variation of site

}

parameters {
 
  //Panicles
  //fixed effects
  real b0_p;    
  real bsize_p;   
  real bsex_p;   
  real bppt_p;  
  real btemp_p;  
  real bcvppt_p;  
  real bppttemp_p;
  real bcvppttemp_p;
  real bsizesex_p;
  real bpptsex_p;
  real btempsex_p;
  real bcvpptsex_p;
  real bppttempsex_p;
  real bcvppttempsex_p;
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
                bsize_p * size_p[ipan] + bsex_p * male_p[ipan] + bppt_p * ppt_p[ipan] + btemp_p * temp_p[ipan]+ bcvppt_p * cvppt_p[ipan] +
                //2-way interactions
                bpptsex_p * ppt_p[ipan] * male_p[ipan] +
                btempsex_p * temp_p[ipan] * male_p[ipan] +
                bcvpptsex_p * cvppt_p[ipan] * male_p[ipan] +
                bsizesex_p * size_p[ipan] * male_p[ipan] +
                bppttemp_p * temp_p[ipan] * ppt_p[ipan]  +
                bcvppttemp_p * cvppt_p[ipan] * temp_p[ipan] +


                //3-way interaction
                bppttempsex_p * temp_p[ipan] * ppt_p[ipan] * male_p[ipan] +
                bcvppttempsex_p * cvppt_p[ipan] * temp_p[ipan]* male_p[ipan] +


                //random effects
                block_rfx_p[block_p[ipan]] +
                source_rfx_p[source_p[ipan]]+
                site_rfx_p[site_p[ipan]];
    }

}

model {
 

  // fertility
  b0_p ~ normal(0, 500);    
  bsize_p ~ normal(0, 100);   
  bsex_p ~ normal(0, 100);   
  bppt_p ~ normal(0, 100);  
  btemp_p ~ normal(0, 100);  
  bcvppt_p ~ normal(0, 100);  
  bpptsex_p ~ normal(0, 100);
  bsizesex_p ~ normal(0, 100);
  bppttemp_p ~ normal(0, 100);
  bcvppttemp_p ~ normal(0, 100);
  bppttempsex_p ~ normal(0, 100); 
  bcvppttempsex_p ~ normal(0, 100);
  btempsex_p ~ normal(0, 100);
  bppttempsex_p ~ normal(0, 100);
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

