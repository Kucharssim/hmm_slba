  // Backward algorithm: log(beta[t, k]) = log p(state_t = k | x_{t:T})
  vector[] backward(int K, int T, vector log_initial, matrix log_transition, vector[] emission_log_likelihoods) {
    vector[K] log_beta[T];

    log_beta[T] = rep_vector(1, K);

    for(tback in 1:(T-1)){
      int t = T-tback;

      for(j in 1:K){
        vector[K] log_current = log_beta[t+1] + log_transition[:,j] + emission_log_likelihoods[t+1, j];
        log_beta[t, j] = log_sum_exp(log_current);
      }
    }

    return log_beta;
  }
