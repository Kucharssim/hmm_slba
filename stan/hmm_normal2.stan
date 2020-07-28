data{
  int<lower=1> N_obs;
  int<lower=1> N_sts;
  real y[N_obs];
}
// transformed data{
//   matrix[N_sts,N_sts] state_tran_prob;
//   real stay_prob = 1.0 - 1.0 / N_sts;
//   
//   for(i in 1:N_sts) {
//     for(j in 1:N_sts){
//       if(i == j) {
//         state_tran_prob[i, j] = stay_prob;
//       } else {
//         state_tran_prob[i, j] = (1.0 - stay_prob) / (N_sts - 1.0);
//       }
//     }
//   }
// }
parameters{
  simplex[N_sts] state_init_prob;
  simplex[N_sts] state_tran_prob[N_sts];
  ordered[N_sts] mu;
  vector<lower=0>[N_sts] sigma;
}
transformed parameters{
  matrix[N_sts, N_obs] emission_log_likelihoods;
  matrix[N_sts, N_sts] state_tran_prob_mat;
  
  for(state in 1:N_sts){
    state_tran_prob_mat[state,] = to_row_vector(state_tran_prob[state]);
    for(n in 1:N_obs){
      emission_log_likelihoods[state, n] = normal_lpdf(y[n] | mu[state], sigma[state]);
    }
    
  }
}
model{
  target += hmm_marginal(emission_log_likelihoods, state_tran_prob_mat, state_init_prob);
  
  state_init_prob ~ dirichlet(rep_vector(2.0, N_sts));
  for(state in 1:N_sts){
    state_tran_prob[state] ~ dirichlet(rep_vector(2, N_sts));
    mu[state]    ~ normal(0, 10);
    sigma[state] ~ gamma(2, 2);
  }
}
generated quantities{
  matrix[N_sts, N_obs] state_prob = hmm_hidden_state_prob(emission_log_likelihoods, state_tran_prob_mat, state_init_prob); 
}

