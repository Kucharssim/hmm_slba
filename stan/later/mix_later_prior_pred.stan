functions{
#include stan/helpers/positive_normal.stan
#include stan/helpers/max.stan
#include stan/later/later_functions.stan
}
data {
  int<lower=0> N_obs; // number of observations
  int<lower=2> N_acc; // number of accumulators
  int<lower=2> N_sts; // number of latent states
  int<lower=1,upper=N_acc> responses[N_obs]; // responses
  real<lower=0> rt[N_obs]; // reaction times

  // copy parameters equated across accumulators
  int<lower=1,upper=N_acc*N_sts> nu_indices   [N_sts,N_acc];
  int<lower=1,upper=N_acc*N_sts> sigma_indices[N_sts,N_acc];
  int<lower=1,upper=N_acc*N_sts> alpha_indices[N_sts,N_acc];
  
  // hyperparameters
  vector<lower=0>[max_int(nu_indices)] nu_alpha;
  
  vector[max_int(sigma_indices)] sigma_mu;
  vector<lower=0>[max_int(sigma_indices)] sigma_sigma;
  
  vector[max_int(alpha_indices)] alpha_mu;
  vector<lower=0>[max_int(alpha_indices)] alpha_sigma;
  
  vector<lower=0>[N_sts] weights_alpha;
}
transformed data{
  int nu_num     = max_int(nu_indices   );
  int sigma_num  = max_int(sigma_indices);
  int alpha_num  = max_int(alpha_indices);
}
parameters {
  simplex[nu_num] nu; // drift rates with sum(nu) = 1 constraint
  vector<lower=0>[sigma_num] sigma; // standard deviations of drift rates
  vector<lower=0>[alpha_num] alpha; // decision boundary
  simplex[N_sts] weights;
}
transformed parameters{
  vector[N_acc] nu_vec   [N_sts];
  vector[N_acc] sigma_vec[N_sts];
  vector[N_acc] alpha_vec[N_sts];
  
  for(k in 1:N_sts){
    nu_vec[k]    = nu[nu_indices[k]];
    sigma_vec[k] = sigma[sigma_indices[k]];
    alpha_vec[k] = alpha[alpha_indices[k]];
    
    nu_vec[k] = nu_vec[k] / sum(nu_vec[k]); // normalize (assure sum == 1 constraint)
  }
}
model{
  for (n in 1:N_obs){
    vector[N_sts] log_lik = log(weights);
    
    for(k in 1:N_sts) log_lik[k] += later_lpdf(rt[n] | responses[n], nu_vec[k], sigma_vec[k], alpha_vec[k], 0);
    
    target += log_sum_exp(log_lik);
  }

  if(include_prior){
    nu ~ dirichlet(nu_alpha);
    for(i in 1:sigma_num) sigma[i] ~ normal(sigma_mu[i], sigma_sigma[i]);
    for(i in 1:alpha_num) alpha[i] ~ normal(alpha_mu[i], alpha_sigma[i]);
    weights ~ dirichlet(weights_alpha);
  }
}
generated quantities{
  int<lower=1, upper=N_sts> state;
  // define data
  real<lower=0> rt[N_obs];
  int<lower=1,upper=N_acc> responses[N_obs];

  // define parameters
  simplex[nu_num] nu = dirichlet_rng(nu_alpha); // drift rates with sum(nu) = 1 constraint
  vector<lower=0>[sigma_num] sigma; // standard deviations of drift rates
  vector<lower=0>[alpha_num] alpha; // decision boundary
  
  // reformat parameters into vectors
  vector[N_acc] nu_vec[N_sts];
  vector[N_acc] sigma_vec[N_sts];
  vector[N_acc] alpha_vec[N_sts];
  
  // draw parameters
  for(i in 1:sigma_num) sigma[i] = positive_normal_rng(sigma_mu[i], sigma_sigma[i]);
  for(i in 1:alpha_num) alpha[i] = positive_normal_rng(alpha_mu[i], alpha_sigma[i]);
  
  fir(k in 1:N_sts){
    nu_vec[k]     = nu   [nu_indices   [k]];
    sigma_vec[k]  = sigma[sigma_indices[k]];
    alpha_vec[k]  = alpha[alpha_indices[k]];
    
    nu_vec[k] = nu_vec[k] / sum(nu_vec[k]); // normalize (assure sum == 1 constraint)
  }
  for(n in 1:N_obs){
    vector[2] pred;
    
    pred = later_rng(nu_vec, sigma_vec, alpha_vec, 0);
    rt[n] = pred[1];
    responses[n] = response_to_integer(pred[2], N_acc);
  }
}