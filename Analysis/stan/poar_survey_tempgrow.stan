

// The input data is a vector 'y' of length 'N'.
data {
  // Data for seed viability sub-model (v)
  int<lower=0> N;   // data points
  int<lower=0> y[N];  // number of female inflorescences
  int<lower=0> n_trials[N]; // number of total inflorescences
  vector[N] tempgrow ; // longitude
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real b0;
  real b_tempgrow;
}

transformed parameters{
  real pred[N];
  for(i in 1:N){
    pred[i] = b0 + b_tempgrow * tempgrow[i];
  }
}

model {
  b0 ~ normal(0, 100);   
  b_tempgrow ~ normal(0, 100);   
  y ~ binomial_logit(n_trials, pred);
}

