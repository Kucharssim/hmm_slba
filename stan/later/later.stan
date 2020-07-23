functions{
#include stan/helpers/positive_normal.stan
#include stan/later/later_functions.stan
}
data {
  int<lower=0> N_obs; // number of observations
  int<lower=2> N_acc; // number of accumulators
  int<lower=1,upper=N_acc> responses[N_obs]; // responses
  real<lower=0> rt[N_obs]; // reaction times

  // copy parameters equated across accumulators
  int<lower=1,upper=N_acc> nu_indices   [N_acc];
  int<lower=1,upper=N_acc> sigma_indices[N_acc];
  int<lower=1,upper=N_acc> alpha_indices[N_acc];
  
  int<lower=0,upper=1> include_prior;
  
  // hyperparameters
  vector<lower=0>[max(nu_indices)] nu_alpha;
  
  vector[max(sigma_indices)] sigma_mu;
  vector<lower=0>[max(sigma_indices)] sigma_sigma;
  
  vector[max(alpha_indices)] alpha_mu;
  vector<lower=0>[max(alpha_indices)] alpha_sigma;
}
transformed data{
  int nu_num     = max(nu_indices   );
  int sigma_num  = max(sigma_indices);
  int alpha_num  = max(alpha_indices);
}
parameters {
  simplex[nu_num] nu; // drift rates with sum(nu) = 1 constraint
  vector<lower=0>[sigma_num] sigma; // standard deviations of drift rates
  vector<lower=0>[alpha_num] alpha; // decision boundary
}
transformed parameters{
  vector[N_acc] nu_vec     = nu   [nu_indices   ];
  vector[N_acc] sigma_vec  = sigma[sigma_indices];
  vector[N_acc] alpha_vec  = alpha[alpha_indices];
  
  nu_vec = nu_vec / sum(nu_vec); // normalize (assure sum == 1 constraint)
}
model{
  for (n in 1:N_obs){
    target += later_lpdf(rt[n] | responses[n], nu_vec, sigma_vec, alpha_vec, 0);
  }

  if(include_prior){
    nu ~ dirichlet(nu_alpha);
    for(i in 1:sigma_num) sigma[i] ~ normal(sigma_mu[i], sigma_sigma[i]);
    for(i in 1:alpha_num) alpha[i] ~ normal(alpha_mu[i], alpha_sigma[i]);
  }
}
