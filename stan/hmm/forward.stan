  // Forward algorithm: log(alpha[t, k]) = log p(state_t = k | x_{1:t})
  vector[] forward(int K, int T, vector log_initial, matrix log_transition, vector[] emission_log_likelihoods) {
    vector[K] log_alpha[T];

    log_alpha[1] = log_initial + emission_log_likelihoods[1];

    for (t in 2:T) {
      for (j in 1:K) {
        vector[K] log_current = log_alpha[t-1] + log_transition[:,j] + emission_log_likelihoods[t, j];
        log_alpha[t, j] = log_sum_exp(log_current);
      }
    }

    return log_alpha;
  }
