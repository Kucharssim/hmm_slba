functions{
#include stan/helpers/positive_normal.stan
#include stan/later/later_functions.stan
}
data {
  int<lower=0> N_obs; // number of observations
  int<lower=2> N_acc; // number of accumulators

  // copy parameters equated across accumulators
  int<lower=1,upper=N_acc> nu_indices   [N_acc];
  int<lower=1,upper=N_acc> sigma_indices[N_acc];
  int<lower=1,upper=N_acc> alpha_indices[N_acc];
  
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
}
model{
}
generated quantities{
  // define data
  real<lower=0> rt[N_obs];
  int<lower=1,upper=N_acc> responses[N_obs];

  // define parameters
  simplex[nu_num] nu = dirichlet_rng(nu_alpha); // drift rates with sum(nu) = 1 constraint
  vector<lower=0>[sigma_num] sigma; // standard deviations of drift rates
  vector<lower=0>[alpha_num] alpha; // decision boundary
  
  // reformat parameters into vectors
  vector[N_acc] nu_vec;
  vector[N_acc] sigma_vec;
  vector[N_acc] alpha_vec;
  
  // draw parameters
  for(i in 1:sigma_num) sigma[i] = positive_normal_rng(sigma_mu[i], sigma_sigma[i]);
  for(i in 1:alpha_num) alpha[i] = positive_normal_rng(alpha_mu[i], alpha_sigma[i]);
  
  nu_vec     = nu   [nu_indices   ];
  sigma_vec  = sigma[sigma_indices];
  alpha_vec  = alpha[alpha_indices];
  
  nu_vec = nu_vec / sum(nu_vec); // normalize (assure sum == 1 constraint)
  
  for(n in 1:N_obs){
    vector[2] pred;
    
    pred = later_rng(nu_vec, sigma_vec, alpha_vec, 0);
    rt[n] = pred[1];
    responses[n] = response_to_integer(pred[2], N_acc);
  }
}
