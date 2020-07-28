functions{
#include stan/helpers/positive_normal.stan
#include stan/helpers/max.stan
#include stan/later/later_functions.stan
}
data {
  int<lower=0> N_obs; // number of observations
  int<lower=2> N_acc; // number of accumulators
  int<lower=2> N_sts; // number of latent states

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
model{
}
generated quantities{
  int<lower=1, upper=N_sts> state[N_obs];
  // define data
  real<lower=0> rt[N_obs];
  int<lower=1,upper=N_acc> responses[N_obs];
  
  // define summary statistics of data
  vector[N_sts] prop_state;
  vector[N_acc] prop_responses[N_sts];
  vector[N_acc] mean_rt       [N_sts];
  // vector[N_acc]   sd_rt       [N_sts];

  // define mixture parameters
  simplex[N_sts] init_prob;
  simplex[N_sts] tran_prob[N_sts];
  
  // define parameters
  simplex[nu_num] nu = dirichlet_rng(nu_alpha); // drift rates with sum(nu) = 1 constraint
  vector<lower=0>[sigma_num] sigma; // standard deviations of drift rates
  vector<lower=0>[alpha_num] alpha; // decision boundary
  vector<lower=0>[t0_num] t0; // non-decision time
  
  // define parameters as arrays of vectors (array[state][accumulator])
  vector[N_acc] nu_vec[N_sts];
  vector[N_acc] sigma_vec[N_sts];
  vector[N_acc] alpha_vec[N_sts];
  real t0_vec[N_sts]; 
  
  // draw parameters
  init_prob = dirichlet_rng(init_prob_alpha);
  for(i in 1:N_sts) tran_prob[i] = dirichlet_rng(tran_prob_alpha[i]);
  
  for(i in 1:sigma_num) sigma[i] = positive_normal_rng(sigma_mu[i], sigma_sigma[i]);
  for(i in 1:alpha_num) alpha[i] = positive_normal_rng(alpha_mu[i], alpha_sigma[i]);
  for(i in 1:t0_num)    t0[i]    = exponential_rng(t0_beta[i]);
  
  for(k in 1:N_sts){
    nu_vec[k]     = nu   [nu_indices   [k]];
    sigma_vec[k]  = sigma[sigma_indices[k]];
    alpha_vec[k]  = alpha[alpha_indices[k]];
    t0_vec[k]     = t0[t0_indices[k]];
    
    nu_vec[k] = nu_vec[k] / sum(nu_vec[k]); // normalize (assure sum == 1 constraint)
    
    prop_state[k]     = 0.0;
    prop_responses[k] = rep_vector(0.0, N_acc);
    mean_rt[k]        = rep_vector(0.0, N_acc);
  }
  
  for(n in 1:N_obs){
    vector[2] pred;
    
    if(n == 1) {
      state[n] = categorical_rng(init_prob);
    } else{
      state[n] = categorical_rng(tran_prob[state[n-1]]);
    }
    
    pred = later_rng(nu_vec[state[n]], sigma_vec[state[n]], alpha_vec[state[n]], t0_vec[state[n]]);
    rt[n] = pred[1];
    responses[n] = response_to_integer(pred[2], N_acc);
    
    // calculate summary stats
    prop_state[state[n]] += 1.0;
    prop_responses[state[n]][responses[n]] += 1.0;
    mean_rt[state[n]][responses[n]] += rt[n];
  }
  
  for(k in 1:N_sts){
    if(prop_state[k] != 0){
      prop_responses[k] /= prop_state[k];
      mean_rt[k] /= prop_state[k];
    }
  }
  prop_state /= N_obs;
}
