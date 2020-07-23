data{
  int<lower=1> N_obs;
  int<lower=1> N_sts;

  vector[N_sts] mu;
  vector<lower=0>[N_sts] sigma;
  
  simplex[N_sts] state_init_prob;
  simplex[N_sts] state_tran_prob[N_sts];
}
parameters{
}
model{
}
generated quantities{
  int  state[N_obs];
  real y[N_obs];

  for(n in 1:N_obs){
    if(n == 1) {
      state[n] = categorical_rng(state_init_prob);
    } else {
      state[n] = categorical_rng(state_tran_prob[state[n-1]]);
    }
    
    y[n] = normal_rng(mu[state[n]], sigma[state[n]]);
  }

}
