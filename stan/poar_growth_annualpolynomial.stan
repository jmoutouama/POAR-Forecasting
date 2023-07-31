// Inverse Gaussian log probability
functions {
  
  real ig_lpdf (real x, real mu, real lambda){
    //vector [num_elements (x)] prob;
    real prob;
    real lprob;
    prob = (lambda/(2*pi()*(x^3)))^0.5*exp(-lambda*(x - mu)^2/(2*mu^2*x));
    //for (i in 1:num_elements(x)) {
      //   prob[i] = (lambda/(2*pi()*(x[i]^3)))^0.5*exp(-lambda*(x[i] - mu)^2/(2*mu^2*x[i]));
      //}
    lprob = log( prob ); 
    //lprob = sum (log(prob)); 
    return lprob;
  }
  
}

data {
  // Data for all vital rates
  int<lower=0> n_sites;         // N. of sites
  int<lower=0> n_sources;         // N. of source populations

  
  // Data for growth sub-model (g)
  int<lower=0> n_g;    // N. of data points for the surival  model
  int<lower=0> n_blocks_g;         // N. of blocks
  int<lower=0> site_g[n_g];  // site index
  int<lower=0> block_g[n_g];  // block index
  int<lower=0> source_g[n_g];  // source index
  int<lower=1> y_g[n_g]; // # tillers at time t+1.
  vector[n_g] size_g;  // log size at time t
  vector[n_g] male_g;  // sex (male=1, fem=0)
  vector[n_g] ppt_g;  // precipitation of site
  vector[n_g] temp_g;  //  temperature of site
  vector[n_g] cvppt_g;  // precipitation coefficient of variation of site

}

parameters {
   
  //Growth
  //fixed effects
  real b0_g;    
  real bsize_g;   
  real bsex_g;  
  real bsizesex_g;  
  real bppt_g;  
  real btemp_g;  
  real bcvppt_g;  
  real bcvppttemp_g;
  real bpptsex_g;
  real btempsex_g;
  real bcvpptsex_g;
  real bppttemp_g;
  real bppttempsex_g;
  real bcvppttempsex_g;
  real bppt2_g;  
  real btemp2_g;
  real bcvppt2_g;
  real bppt2sex_g;
  real btemp2sex_g;
  real bcvppt2sex_g;
  //random effects
  real<lower=0> block_tau_g; 
  real block_rfx_g[n_blocks_g];  
  real<lower=0> source_tau_g; 
  real source_rfx_g[n_sources];
  real<lower=0> site_tau_g; 
  real site_rfx_g[n_sites];
  real<lower=0> sigma;      // IG shape
  real<lower=0> theta[n_g]; //observation-level deviates

  }

transformed parameters {
  real predG[n_g];
  // prediction for  growth
  for(igrow in 1:n_g){
    predG[igrow] = exp(b0_g + 
                //main effects
                bsize_g * size_g[igrow] + bsex_g * male_g[igrow] + bppt_g * ppt_g[igrow] + btemp_g * temp_g[igrow]+ bcvppt_g * cvppt_g[igrow] +
                //2-way interactions
                bpptsex_g * ppt_g[igrow] * male_g[igrow] +
                btempsex_g * temp_g[igrow] * male_g[igrow] +
                bcvpptsex_g * cvppt_g[igrow] * male_g[igrow] +
                bcvppttemp_g * cvppt_g[igrow] * temp_g[igrow] +
                bsizesex_g * size_g[igrow] * male_g[igrow] +

                //3-way interaction
                bppttempsex_g * temp_g[igrow] * ppt_g[igrow] * male_g[igrow] +
                bcvppttempsex_g * cvppt_g[igrow] * temp_g[igrow]* male_g[igrow] +

                //polynomial 2

                bppt2_g * pow(ppt_g [igrow],2) + 
                btemp2_g * pow(temp_g [igrow],2) + 
                bcvppt2_g * pow(cvppt_g [igrow],2) + 
                bppt2sex_g * male_g [igrow] * pow(ppt_g [igrow],2) +
                btemp2sex_g * male_g [igrow] * pow(temp_g[igrow],2) +
                bcvppt2sex_g * male_g [igrow] * pow(cvppt_g[igrow],2) +
               
                

                //random effects
                block_rfx_g[block_g[igrow]] +
                source_rfx_g[source_g[igrow]]+
                site_rfx_g[site_g[igrow]]);
    }
  
}

model {
  // priors on parameters
  // Growth
  b0_g ~ normal(0, 100);    
  bsize_g ~ normal(0, 100);   
  bsex_g ~ normal(0, 100);   
  bppt_g ~ normal(0, 100);  
  btemp_g ~ normal(0, 100);  
  bcvppt_g ~ normal(0, 100); 
  bppttemp_g ~ normal(0, 100); 
  bcvppttemp_g ~ normal(0, 100);
  bpptsex_g ~ normal(0, 100);
  bsizesex_g ~ normal(0, 100);
  btempsex_g ~ normal(0, 100);
  bppttempsex_g ~ normal(0, 100);
  bcvppttempsex_g ~ normal(0, 100);
  bppt2_g ~ normal(0, 1);
  btemp2_g ~ normal(0, 1);
  bcvppt2_g ~ normal(0, 1);
  bppt2sex_g ~ normal(0, 1);
  btemp2sex_g ~ normal(0, 1);
  bcvppt2sex_g ~ normal(0, 1);
  block_tau_g ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_blocks_g){
    block_rfx_g[i] ~ normal(0, block_tau_g);
  }
  source_tau_g ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_sources){
    source_rfx_g[i] ~ normal(0, source_tau_g);
  }
  site_tau_g ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_sites){
    site_rfx_g[i] ~ normal(0, site_tau_g);
  }
  for(i in 1:n_g){
    theta[i] ~ ig(1, sigma);
  }
  

  // sampling  

  //growth needs loop for zero truncation
  for (i in 1:n_g) {
    //T[1,] for zero truncation
    y_g[i] ~ poisson(predG[i] * theta[i]) T[1,]; 
  }
  
}

generated quantities {
  vector[n_g] log_lik;
  for (i in 1:n_g) {
  log_lik[i] = poisson_lpmf (y_g[i] | predG[i]);
 }
}


