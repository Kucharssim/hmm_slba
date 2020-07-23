functions{
#include hmm/forward.stan
#include hmm/backward.stan
#include hmm/forwardbackward.stan
}
data{
  int<lower=1> N_obs;
  int<lower=1> N_sts;
  real y[N_obs];
}
parameters{
  ordered[N_sts] mu;
  vector<lower=0>[N_sts] sigma;
  
  simplex[N_sts] state_init_prob;
  simplex[N_sts] state_tran_prob[N_sts];
}
transformed parameters{
  vector[N_sts] log_state_init_prob = log(state_init_prob);
  matrix[N_sts,N_sts] log_state_tran_prob;
  vector[N_sts] emission_log_likelihoods[N_obs];
  vector[N_sts] log_alpha[N_obs]; // forward variable
  
  for(state in 1:N_sts){
    log_state_tran_prob[state, ] = to_row_vector(log(state_tran_prob[state]));
  }
  
  for(n in 1:N_obs){
    for(state in 1:N_sts){
      emission_log_likelihoods[n, state] = normal_lpdf(y[n] | mu[state], sigma[state]);
    }
  }

  log_alpha = forward(N_sts, N_obs, log_state_init_prob, log_state_tran_prob, emission_log_likelihoods);
}
model{
  target += log_sum_exp(log_alpha[N_obs]);
  
  state_init_prob ~ dirichlet(rep_vector(2, N_sts));
  for(state in 1:N_sts){
    state_tran_prob[state] ~ dirichlet(rep_vector(2, N_sts));
    mu[state]              ~ normal(0, 10);
    sigma[state]           ~ gamma(2, 2);
  }
}
generated quantities{
  vector [N_sts] log_beta[N_obs]; // backward variable
  simplex[N_sts] gamma[N_obs]; // gamma[t, k] = p(state_t = k | x_{1:T})
  real y_rep[N_obs]; // posterior predictives

  log_beta = backward       (N_sts, N_obs, log_state_init_prob, log_state_tran_prob, emission_log_likelihoods);
  gamma    = forwardbackward(N_sts, N_obs, log_state_init_prob, log_state_tran_prob, emission_log_likelihoods, log_alpha, log_beta);

  for(n in 1:N_obs){
    int state = categorical_rng(gamma[n]);
    
    y_rep[n] = normal_rng(mu[state], sigma[state]);
  }

}

