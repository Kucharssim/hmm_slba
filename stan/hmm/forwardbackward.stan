  // Forward-backward algorithm: gamma[t, k] = p(state_t = k | x_{1:T})
  vector[] forwardbackward(int K, int T, vector log_initial, matrix log_transition, vector[] emission_log_likelihoods, vector[] log_alpha, vector[] log_beta) {
    vector[K] log_gamma[T];
    vector[K] gamma[T];

    for(t in 1:T){
      log_gamma[t] = log_alpha[t] + log_beta[t];
      gamma[t] = softmax(log_gamma[t]);
    }

    return gamma;
  } // Forward-backward
