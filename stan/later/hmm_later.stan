functions{
#include stan/helpers/positive_normal.stan
#include stan/helpers/max.stan
#include stan/helpers/array_to_matrix.stan
#include stan/later/later_functions.stan
}
data {
  int<lower=0> N_obs; // number of observations
  int<lower=2> N_acc; // number of accumulators
  int<lower=2> N_sts; // number of latent states
  real<lower=0> rt[N_obs];
  int<lower=1,upper=N_acc> responses[N_obs];

  // copy parameters equated across accumulators
  int<lower=1,upper=N_acc*N_sts> nu_indices   [N_sts,N_acc];
  int<lower=1,upper=N_acc*N_sts> sigma_indices[N_sts,N_acc];
  int<lower=1,upper=N_acc*N_sts> alpha_indices[N_sts,N_acc];
  int<lower=1,upper=N_sts>       t0_indices[N_sts];
  
  // hyperparameters
  vector<lower=0>[max_int(nu_indices)] nu_alpha;
  
  vector[max_int(sigma_indices)] sigma_mu;
  vector<lower=0>[max_int(sigma_indices)] sigma_sigma;
  
  vector[max_int(alpha_indices)] alpha_mu;
  vector<lower=0>[max_int(alpha_indices)] alpha_sigma;
  
  vector<lower=0>[max(t0_indices)] t0_beta;
  
  vector<lower=0>[N_sts] init_prob_alpha;
  vector<lower=0>[N_sts] tran_prob_alpha[N_sts];
  
}
transformed data{
  int nu_num     = max_int(nu_indices   );
  int sigma_num  = max_int(sigma_indices);
  int alpha_num  = max_int(alpha_indices);
  int t0_num     = max(t0_indices);
}
parameters{
  // define parameters
  simplex[N_sts] init_prob;
  simplex[N_sts] tran_prob[N_sts];
  
  simplex[nu_num] nu; // drift rates with sum(nu) = 1 constraint
  vector<lower=0>[sigma_num] sigma; // standard deviations of drift rates
  vector<lower=0>[alpha_num] alpha; // decision boundary
  vector<lower=0>[t0_num] t0; // non-decision time
}
transformed parameters{
  // define parameters as arrays of vectors (array[state][accumulator])
  vector[N_acc] nu_vec[N_sts];
  vector[N_acc] sigma_vec[N_sts];
  vector[N_acc] alpha_vec[N_sts];
  real t0_vec[N_sts]; 
  
  matrix[N_sts, N_obs] emission_log_likelihoods; // log-likelihoods of observations under states
  
  for(k in 1:N_sts){
    nu_vec[k]     = nu   [nu_indices   [k]];
    sigma_vec[k]  = sigma[sigma_indices[k]];
    alpha_vec[k]  = alpha[alpha_indices[k]];
    t0_vec[k]     = t0[t0_indices[k]];
    
    nu_vec[k] = nu_vec[k] / sum(nu_vec[k]); // normalize (assure sum == 1 constraint)
  }
  
  for(k in 1:N_sts){
    for(n in 1:N_obs){
      emission_log_likelihoods[k, n] = later_lpdf(rt[n] | responses[n], nu_vec[k], sigma_vec[k], alpha_vec[k], t0_vec[k]);
    }
  }
  
}
model{
  // likelihood
  target += hmm_marginal(emission_log_likelihoods, array_to_matrix(tran_prob, N_sts), init_prob);
  
  // priors
  init_prob ~ dirichlet(init_prob_alpha);
  for(i in 1:N_sts) tran_prob[i] ~ dirichlet(tran_prob_alpha[i]);
  
  nu ~ dirichlet(nu_alpha);
  for(i in 1:sigma_num) sigma[i] ~ normal(sigma_mu[i], sigma_sigma[i]);
  for(i in 1:alpha_num) alpha[i] ~ normal(alpha_mu[i], alpha_sigma[i]);
  for(i in 1:t0_num)    t0[i]    ~ exponential(t0_beta[i]);
}
generated quantities{
  matrix[N_sts, N_obs] state_prob = hmm_hidden_state_prob(emission_log_likelihoods, array_to_matrix(tran_prob, N_sts), init_prob); 
}
