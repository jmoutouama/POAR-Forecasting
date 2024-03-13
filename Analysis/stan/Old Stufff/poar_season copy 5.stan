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

  // Data for growth sub-model (g)
  int<lower=0> n_g;    // N. of data points for the surival  model
  int<lower=0> n_blocks_g;         // N. of blocks
  int<lower=0> site_g[n_g];  // site index
  int<lower=0> block_g[n_g];  // block index
  int<lower=0> source_g[n_g];  // source index
  int<lower=1> y_g[n_g]; // number tillers at time t+1.
  vector[n_g] size_g;  //   size at time t
  vector[n_g] male_g;  // sex (male=1, fem=0)
  vector[n_g] pptgrow_g;  // precipitation of growing season 
  vector[n_g] pptdorm_g;  // precipitation of dormant season
  vector[n_g] tempgrow_g;  //  temperature of growing season
  vector[n_g] tempdorm_g;  //  temperature of growing season

  // Data for flowering sub-model (f)
  int<lower=0> n_f;    // N. of data points for the flowering  model
  int<lower=0> n_blocks_f;         // N. of blocks
  int<lower=0> site_f[n_f];  // site index
  int<lower=0> block_f[n_f];  // block index
  int<lower=0> source_f[n_f];  // source index
  int<lower=0,upper=1> y_f[n_f]; // Flowering at time t+1.
  vector[n_f] size_f;  // log size at time t
  vector[n_f] male_f;  // sex (male=1, fem=0)
  vector[n_f] pptgrow_f;  // precipitation of growing season 
  vector[n_f] pptdorm_f;  // precipitation of dormant season
  vector[n_f] tempgrow_f;  //  temperature of growing season
  vector[n_f] tempdorm_f;  //  temperature of growing season

  // Data for fertility sub-model (p)
  int<lower=0> n_p;    // N. of data points for the fertility model
  int<lower=0> n_blocks_p;         // N. of blocks
  int<lower=0> site_p[n_p];  // site index
  int<lower=0> block_p[n_p];  // block index
  int<lower=0> source_p[n_p];  // source index
  int<lower=1> y_p[n_p]; // number panicles at time t.
  vector[n_p] size_p;  // log size at time t
  vector[n_p] male_p;  // sex (male=1, fem=0)
  vector[n_p] pptgrow_p;  // precipitation of growing season 
  vector[n_p] pptdorm_p;  // precipitation of dormant season
  vector[n_p] tempgrow_p;  //  temperature of growing season
  vector[n_p] tempdorm_p;  //  temperature of growing season

  // Data for seed viability sub-model (v)
  int<lower=0> n_v;   // data points
  int<lower=0> y_v[n_v];  // number of viable seeds
  int<lower=0> tot_seeds_v[n_v]; // number of trials
  real SR_v[n_v]; // Sex ratio (proportion female)
  
  // Data for seed germination sub-model (m)
  int<lower=0> n_m;   // data points
  int<lower=0> y_m[n_m];  // number of germinating seeds
  int<lower=0> tot_seeds_m[n_m]; // number of trials
  real SR_m[n_m]; // Sex ratio (proportion female)
  
  // data for seed count
  int<lower=0> n_d;   // data points
  int<lower=0> y_d[n_d];  // total number of seeds per panicle 
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
  real bsizesex_s;
  real btempdormpptdorm_s;
  real btempgrowpptgrow_s;
  real bpptgrowsex_s;  
  real bpptdormsex_s; 
  real btempgrowsex_s;  
  real btempdormsex_s;  
  real btempdormpptdormsex_s;  
  real btempgrowpptgrowsex_s;  
  real bpptgrow2_s;  
  real bpptdorm2_s;  
  real btempgrow2_s;  
  real btempdorm2_s; 
  real bpptgrow2sex_s;  
  real bpptdorm2sex_s; 
  real btempgrow2sex_s;  
  real btempdorm2sex_s; 
  //random effects
  real<lower=0> block_tau_s; 
  real block_rfx_s[n_blocks_s];  
  real<lower=0> source_tau_s; 
  real source_rfx_s[n_sources];
  real<lower=0> site_tau_s; 
  real site_rfx_s[n_sites];

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

  //Flowering
  //fixed effects
  real b0_f;    
  real bsize_f;   
  real bsizesex_f;   
  real bsex_f;   
  real bpptgrow_f;  
  real bpptdorm_f;  
  real btempgrow_f;  
  real btempdorm_f;  
  real bpptgrowsex_f;  
  real bpptdormsex_f; 
  real btempgrowsex_f;  
  real btempdormsex_f;  
  real btempdormpptdorm_f;  
  real btempgrowpptgrow_f;  
  real btempdormpptdormsex_f;  
  real btempgrowpptgrowsex_f;  
  real bpptgrow2_f;  
  real bpptdorm2_f;  
  real btempgrow2_f;  
  real btempdorm2_f;  
  real bpptgrow2sex_f;  
  real bpptdorm2sex_f; 
  real btempgrow2sex_f;  
  real btempdorm2sex_f; 
  //random effects
  real<lower=0> block_tau_f; 
  real block_rfx_f[n_blocks_f];  
  real<lower=0> source_tau_f; 
  real source_rfx_f[n_sources];
  real<lower=0> site_tau_f; 
  real site_rfx_f[n_sites];

  //Panicles
  //fixed effects
  real b0_p;    
  real bsize_p;   
  real bsex_p;   
  real bpptgrow_p;  
  real bpptdorm_p;  
  real btempgrow_p;  
  real btempdorm_p;  
  real bsizesex_p;
  real bpptgrowsex_p;  
  real bpptdormsex_p; 
  real btempgrowsex_p;  
  real btempdormsex_p;  
  real btempdormpptdormsex_p;  
  real btempgrowpptgrowsex_p;  
  real btempdormpptdorm_p;
  real btempgrowpptgrow_p;
  real bpptgrow2_p;  
  real bpptdorm2_p;  
  real btempgrow2_p;  
  real btempdorm2_p; 
  real bpptgrow2sex_p;  
  real bpptdorm2sex_p; 
  real btempgrow2sex_p;  
  real btempdorm2sex_p; 
  //random effects
  real<lower=0> block_tau_p; 
  real block_rfx_p[n_blocks_p];  
  real<lower=0> source_tau_p; 
  real source_rfx_p[n_sources];
  real<lower=0> site_tau_p; 
  real site_rfx_p[n_sites];
  real<lower=0> phi_p; // Panicle dispersion parameter

  // Seed viability parameters
  real<lower=0,upper=1> v0;
  real<lower=0> a_v;
  real<lower=0.1> phi_v;             
  // Germination rate
  real<lower=0,upper=1> m;
  real<lower=0.1> phi_m;  
  // Seed count
  real<lower=0> lambda_d;
  }

transformed parameters {

  real predS[n_s];
  real predG[n_g];
  real predF[n_f];
  real predP[n_p];
  real<lower=0,upper=1> predV[n_v];
  real<lower=0,upper=1> predM[n_m];

   // beta-binom reparameterization for viab and germ
  vector[n_v] alpha_v; 
  vector[n_v] beta_v;   
  vector[n_m] alpha_m; 
  vector[n_m] beta_m;  

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
                bsizesex_s * size_s [isurv] * male_s[isurv] +
                btempdormpptdorm_s * tempdorm_s[isurv] * pptdorm_s[isurv] +
                btempgrowpptgrow_s * tempgrow_s[isurv] * pptgrow_s[isurv] +
                //3-way interaction
                btempdormpptdormsex_s * tempdorm_s[isurv] * pptdorm_s[isurv] * male_s[isurv] +
                btempgrowpptgrowsex_s * tempgrow_s[isurv] * tempgrow_s[isurv] * male_s[isurv] +
                //polynomial 2
                bpptgrow2_s * pow(pptgrow_s[isurv],2) + 
                bpptdorm2_s * pow(pptdorm_s[isurv],2) + 
                btempgrow2_s *  pow(tempgrow_s[isurv],2) + 
                btempdorm2_s * pow(tempdorm_s[isurv],2) + 
                bpptgrow2sex_s * male_s[isurv] * pow(pptgrow_s[isurv],2) + 
                bpptdorm2sex_s * male_s[isurv] * pow(pptdorm_s[isurv],2) + 
                btempgrow2sex_s *  male_s[isurv] * pow(tempgrow_s[isurv],2) + 
                btempdorm2sex_s * male_s[isurv] * pow(tempdorm_s[isurv],2) + 
                //random effects
                block_rfx_s[block_s[isurv]] +
                source_rfx_s[source_s[isurv]]+
                site_rfx_s[site_s[isurv]];
    }
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

  // prediction for flowering
  for(iflow in 1:n_f){
    predF[iflow] = b0_f + 
                //main effects
                bsize_f * size_f[iflow] + bsex_f * male_f[iflow] + bpptgrow_f * pptgrow_f[iflow] + bpptdorm_f * pptdorm_f[iflow] + btempgrow_f * tempgrow_f[iflow]  + btempdorm_f * tempdorm_f [iflow] +
                //2-way interactions
                bpptgrowsex_f * pptgrow_f[iflow] * male_f[iflow] +
                bpptdormsex_f * pptdorm_f[iflow] * male_f[iflow] +
                btempgrowsex_f * tempgrow_f[iflow] * male_f[iflow] +
                btempdormsex_f * tempdorm_f[iflow] * male_f[iflow] +
                bsizesex_f * size_f[iflow] * male_f[iflow] +
                btempdormpptdorm_f * tempdorm_f[iflow] * pptdorm_f[iflow]  +
                btempgrowpptgrow_f * tempgrow_f[iflow] * pptgrow_f[iflow]  +
                //3-way interaction
                btempdormpptdormsex_f * tempdorm_f[iflow] * pptdorm_f[iflow] * male_f[iflow] +
                btempgrowpptgrowsex_f * tempgrow_f[iflow] * pptgrow_f[iflow] * male_f[iflow] +
                //polynomial 2
                bpptgrow2_f * pow(pptgrow_f[iflow],2) + 
                bpptdorm2_f * pow(pptdorm_f[iflow],2) + 
                btempgrow2_f * pow(tempgrow_f[iflow],2) + 
                btempdorm2_f * pow(tempdorm_f[iflow],2) + 
                bpptgrow2sex_f * male_f[iflow] * pow(pptgrow_f[iflow],2) + 
                bpptdorm2sex_f * male_f[iflow] * pow(pptdorm_f[iflow],2) + 
                btempgrow2sex_f *  male_f[iflow] * pow(tempgrow_f[iflow],2) + 
                btempdorm2sex_f * male_f[iflow] * pow(tempdorm_f[iflow],2) + 
                //random effects
                block_rfx_f[block_f[iflow]] +
                source_rfx_f[source_f[iflow]]+
                site_rfx_f[site_f[iflow]];
    }

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
                //polynomial 2
                bpptgrow2_p * pow(pptgrow_p[ipan],2) + 
                bpptdorm2_p * pow(pptdorm_p[ipan],2) + 
                btempgrow2_p * pow(tempgrow_p[ipan],2) + 
                btempdorm2_p * pow(tempdorm_p[ipan],2) + 
                bpptgrow2sex_p * male_p[ipan] * pow(pptgrow_p[ipan],2) + 
                bpptdorm2sex_p * male_p[ipan] * pow(pptdorm_p[ipan],2) + 
                btempgrow2sex_p *  male_p[ipan] * pow(tempgrow_p[ipan],2) + 
                btempdorm2sex_p * male_p[ipan] * pow(tempdorm_p[ipan],2) + 
                //random effects
                block_rfx_p[block_p[ipan]] +
                source_rfx_p[source_p[ipan]]+
                site_rfx_p[site_p[ipan]];
    }

  // Prediction for seed viability
  for(iviab in 1:n_v){
    predV[iviab] = v0 * (1 - pow(SR_v[iviab],a_v) ) + 0.00001;
    alpha_v[iviab] = predV[iviab] * phi_v;
    beta_v[iviab] = (1 - predV[iviab]) * phi_v; 
    }
  // Prediction for germination
  for(igerm in 1:n_m){
    predM[igerm] = m * v0 * (1 - pow(SR_m[igerm],a_v) ) + 0.00001;
    alpha_m[igerm] = predM[igerm] * phi_m;
    beta_m[igerm] = (1 - predM[igerm]) * phi_m; 
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
  bsizesex_s ~ normal(0, 100); 
  btempdormpptdorm_s ~ normal(0, 100); 
  btempgrowpptgrow_s ~ normal(0, 100); 
  btempdorm_s ~ normal(0, 100);
  bpptgrowsex_s ~ normal(0, 100);
  bpptdormsex_s ~ normal(0, 100);
  btempgrowsex_s ~ normal(0, 100);
  btempdormsex_s ~ normal(0, 100);
  btempdormpptdormsex_s ~ normal(0, 100);
  btempgrowpptgrowsex_s ~ normal(0, 100);
  bpptgrow2sex_s ~ normal(0, 1);  
  bpptdorm2sex_s ~ normal(0, 1); 
  btempgrow2sex_s ~ normal(0, 1);  
  btempdorm2sex_s ~ normal(0, 1); 
  bpptgrow2_s ~ normal(0, 1); 
  bpptdorm2_s ~ normal(0, 1); 
  btempgrow2_s ~ normal(0, 1); 
  btempdorm2_s ~ normal(0, 1); 
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

  // Flowering
  b0_f ~ normal(0, 500);    
  bsize_f ~ normal(0, 100);  
  bsizesex_f ~ normal(0, 100);  
  bsex_f ~ normal(0, 100);   
  bpptgrow_f ~ normal(0, 100);  
  bpptdorm_f ~ normal(0, 100);  
  btempgrow_f ~ normal(0, 100);  
  btempdorm_f ~ normal(0, 100);
  bpptgrowsex_f ~ normal(0, 100);
  bpptdormsex_f ~ normal(0, 100);
  btempgrowsex_f ~ normal(0, 100);
  btempdormsex_f ~ normal(0, 100);
  btempdormpptdorm_f ~ normal(0, 100);
  btempgrowpptgrow_f ~ normal(0, 100);
  btempdormpptdormsex_f ~ normal(0, 1);
  btempgrowpptgrowsex_f ~ normal(0, 1);
  bpptgrow2_f ~ normal(0, 1); 
  bpptdorm2_f ~ normal(0, 1);
  btempgrow2_f ~ normal(0, 1);
  btempdorm2_f ~ normal(0, 1);
  bpptgrow2sex_f ~ normal(0, 1);
  bpptdorm2sex_f ~ normal(0, 1); 
  btempgrow2sex_f ~ normal(0, 1); 
  btempdorm2sex_f ~ normal(0, 1); 
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
  // fertility
  b0_p ~ normal(0, 500);    
  bsize_p ~ normal(0, 100);   
  bsex_p ~ normal(0, 100);   
  bpptgrow_p ~ normal(0, 100);  
  bpptdorm_p ~ normal(0, 100);  
  btempgrow_p ~ normal(0, 100);  
  btempdorm_p ~ normal(0, 100);
  bsizesex_p ~ normal(0, 100);
  btempdormpptdorm_p ~ normal(0, 100);
  btempgrowpptgrow_p ~ normal(0, 100);
  bpptgrowsex_p ~ normal(0, 100);
  bpptdormsex_p ~ normal(0, 100);
  btempgrowsex_p ~ normal(0, 100);
  btempdormsex_p ~ normal(0, 100);
  btempdormpptdormsex_p ~ normal(0, 100);
  btempgrowpptgrowsex_p ~ normal(0, 100);
  bpptgrow2sex_p ~ normal(0, 0.5);  
  bpptdorm2sex_p ~ normal(0, 0.5); 
  btempgrow2sex_p ~ normal(0, 0.5);  
  btempdorm2sex_p ~ normal(0, 0.5); 
  bpptgrow2_p ~ normal(0, 0.5); 
  bpptdorm2_p ~ normal(0, 0.5); 
  btempgrow2_p ~ normal(0, 0.5); 
  btempdorm2_p ~ normal(0, 0.5); 
  block_tau_p ~ inv_gamma(0.1, 0.1);
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

  // viability and germination
  v0  ~ beta(10,1);  // intercept viability model
  a_v ~ gamma(1,1);  // "decay" in viability with SR
  phi_v ~ pareto(0.1,1.5);
  m  ~ beta(10,1);  // intercept viability model
  phi_m ~ pareto(0.1,1.5);
  lambda_d ~ inv_gamma(0.001, 0.001);

  // sampling  

  //survival
  y_s ~ bernoulli_logit(predS);

  //growth needs loop for zero truncation
  for (i in 1:n_g) {
    //T[1,] for zero truncation
    y_g[i] ~ poisson(predG[i] * theta[i]) T[1,]; 
  }

  //flowering
  y_f ~ bernoulli_logit(predF);

  //fertility need loop for zero truncation
  for (i in 1:n_p) {
    y_p[i] ~ neg_binomial_2_log(predP[i], phi_p);
    // manually zero-truncating
    target += - log1m(neg_binomial_2_log_lpmf(0 | predP[i], phi_p)); 
  }

  //viability
  y_v ~ beta_binomial(tot_seeds_v, alpha_v, beta_v);
  //germination
  y_m ~ beta_binomial(tot_seeds_m, alpha_m, beta_m);
  //seed number
  y_d ~ poisson(lambda_d);
}




