
data {

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
}

parameters {
    // Seed viability parameters
  real<lower=0,upper=1> v0;
  real<lower=0> a_v;
  real<lower=0.1> phi_v;             
  // Germination rate
  real<lower=0,upper=1> m;
  real<lower=0.1> phi_m;  
}

transformed parameters {

  real<lower=0,upper=1> predV[n_v];
  real<lower=0,upper=1> predM[n_m];
  
  // beta-binom reparameterization for viab and germ
  vector[n_v] alpha_v; 
  vector[n_v] beta_v;   
  vector[n_m] alpha_m; 
  vector[n_m] beta_m;  

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
  // viability and germination
  v0  ~ beta(10,1);  // intercept viability model
  a_v ~ gamma(1,1);  // "decay" in viability with SR
  phi_v ~ pareto(0.1,1.5);
  m  ~ beta(10,1);  // intercept viability model
  phi_m ~ pareto(0.1,1.5);

  // sampling  
  //viability
  y_v ~ beta_binomial(tot_seeds_v, alpha_v, beta_v);
  //germination
  y_m ~ beta_binomial(tot_seeds_m, alpha_m, beta_m);
  
}



