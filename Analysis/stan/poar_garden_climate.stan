

// The input data is a vector 'y' of length 'N'.
data {
  // Data for seed viability sub-model (v)
  int<lower=0> N;   // data points
  int<lower=0> y[N];  // number of female inflorescences
  int<lower=0> n_trials[N]; // number of total inflorescences
  vector[N] pptgrow;  // precipitation of growing season 
  vector[N] pptdorm;  // precipitation of dormant season
  vector[N] tempgrow;  //  temperature of growing season
  vector[N] tempdorm;  //  temperature of growing season
}

// The parameters accepted by the model. Our model
parameters {
  real b0;
  real b_pptgrow;  
  real b_pptdorm;  
  real b_tempgrow;  
  real b_tempdorm; 
  real b_tempdormpptdorm; 
  real b_tempgrowpptgrow; 
  real b_pptgrow2;  
  real b_pptdorm2;  
  real b_tempgrow2;  
  real b_tempdorm2;
}

transformed parameters{
  real pred[N];
  for(i in 1:N){
  // main effects
    pred[i] = b0 + b_pptdorm * pptdorm[i] + b_tempdorm * tempdorm[i] + b_pptgrow * pptgrow[i] + b_tempgrow * tempgrow[i] + 
    // 2-way interactions
    b_tempdormpptdorm * tempdorm[i] * pptdorm[i] +
    b_tempgrowpptgrow * tempgrow[i] * pptgrow[i] +
    //polynomial 2
    b_pptgrow2 * pow(pptgrow[i],2) + 
    b_pptdorm2 * pow(pptdorm[i],2) + 
    b_tempgrow2 *  pow(tempgrow[i],2) + 
    b_tempdorm2 * pow(tempdorm[i],2); 
  }
}

model {
  b0 ~ normal(0, 100);   
  b_pptgrow ~ normal(0,100);  
  b_pptdorm ~ normal (0,100);  
  b_tempgrow ~ normal (0,100);  
  b_tempdorm ~ normal (0,100); 
  b_pptgrow2 ~ normal (0,100);  
  b_pptdorm2 ~ normal (0,100);  
  b_tempgrow2 ~ normal (0,100);  
  b_tempdorm2 ~ normal (0,100);
  b_tempdormpptdorm ~ normal (0,100);
  b_tempgrowpptgrow ~ normal (0,100);
  y ~ binomial_logit(n_trials, pred);
}

