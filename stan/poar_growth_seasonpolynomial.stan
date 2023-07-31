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
  vector[n_g] pptgrow_g;  // precipitation of growing season 
  vector[n_g] pptdorm_g;  // precipitation of dormant season
  vector[n_g] tempgrow_g;  //  temperature of growing season
  vector[n_g] tempdorm_g;  //  temperature of growing season
}

parameters {
   
  //Growth
  //fixed effects
  real b0_g;    
  real bsize_g;   
  real bsizesex_g;   
  real bsex_g;   
  real bpptgrow_g;  
  real bpptdorm_g;  
  real btempgrow_g;  
  real btempdorm_g;  
  real bpptgrowsex_g;  
  real bpptdormsex_g; 
  real btempgrowsex_g;  
  real btempdormsex_g;  
  real btempdormpptdorm_g;  
  real btempgrowpptgrow_g;  
  real btempdormpptdormsex_g;  
  real btempgrowpptgrowsex_g;  
  real bpptgrow2_g;  
  real bpptdorm2_g;  
  real btempgrow2_g;  
  real btempdorm2_g; 
  real bpptgrow2sex_g;  
  real bpptdorm2sex_g; 
  real btempgrow2sex_g;  
  real btempdorm2sex_g; 


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
                bsize_g * size_g[igrow] + bsex_g * male_g[igrow] + bpptgrow_g * pptgrow_g[igrow] + bpptdorm_g * pptdorm_g[igrow] + btempgrow_g * tempgrow_g[igrow]  + btempdorm_g * tempdorm_g [igrow] +
                
                //2-way interactions
                bpptgrowsex_g * pptgrow_g[igrow] * male_g[igrow] +
                bpptdormsex_g * pptdorm_g[igrow] * male_g[igrow] +
                btempgrowsex_g * tempgrow_g[igrow] * male_g[igrow] +
                btempdormsex_g * tempdorm_g[igrow] * male_g[igrow] +
                btempdormpptdorm_g * tempdorm_g[igrow] * pptdorm_g[igrow]  +
                btempgrowpptgrow_g * tempgrow_g[igrow] * pptgrow_g [igrow] +
                bsizesex_g * size_g[igrow] *  male_g[igrow] +

                //3-way interaction
                btempdormpptdormsex_g * tempdorm_g[igrow] * pptdorm_g[igrow] * male_g[igrow] +
                btempgrowpptgrowsex_g * tempgrow_g[igrow] * pptgrow_g[igrow] * male_g[igrow] +

                //polynomial 2
                bpptgrow2_g * pow(pptgrow_g[igrow],2) + 
                bpptdorm2_g * pow(pptdorm_g[igrow],2) + 
                btempgrow2_g * pow(tempgrow_g[igrow],2) + 
                btempdorm2_g * pow(tempdorm_g[igrow],2) + 
                bpptgrow2sex_g * male_g[igrow] * pow(pptgrow_g[igrow],2) + 
                bpptdorm2sex_g * male_g[igrow] * pow(pptdorm_g[igrow],2) + 
                btempgrow2sex_g *  male_g[igrow] * pow(tempgrow_g[igrow],2) + 
                btempdorm2sex_g * male_g[igrow] * pow(tempdorm_g[igrow],2) + 
                
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
  bpptgrow_g ~ normal(0, 100);  
  bpptdorm_g ~ normal(0, 100);  
  btempgrow_g ~ normal(0, 100);  
  btempdorm_g ~ normal(0, 100);
  bsizesex_g ~ normal(0, 100);
  btempdormpptdorm_g ~ normal(0, 100);
  btempgrowpptgrow_g ~ normal(0, 100);
  bpptgrowsex_g ~ normal(0, 100);
  bpptdormsex_g ~ normal(0, 100);
  btempgrowsex_g ~ normal(0, 100);
  btempdormsex_g ~ normal(0, 100);
  btempdormpptdormsex_g ~ normal(0, 100);
  btempgrowpptgrowsex_g ~ normal(0, 100);
  bpptgrow2sex_g ~ normal(0, 0.5);  
  bpptdorm2sex_g ~ normal(0, 0.5); 
  btempgrow2sex_g ~ normal(0, 0.5);  
  btempdorm2sex_g ~ normal(0, 0.5); 
  bpptgrow2_g ~ normal(0, 0.5); 
  bpptdorm2_g ~ normal(0, 0.5); 
  btempgrow2_g ~ normal(0, 0.5); 
  btempdorm2_g ~ normal(0, 0.5); 
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

