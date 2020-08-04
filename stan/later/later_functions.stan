// depends on positive_normal.stan functions
// functions{
  /*
  SINGLE ACCUMULATOR FUNCTIONS
  */
  // rng of a t of a single accumulator
  real later_single_rng(real nu, real sigma, real alpha){
    real y = positive_normal_rng(nu, sigma);
    
    return alpha/y;
  }
  // pdf of a t of a single accumulator
  real later_single_lpdf(real t, real nu, real sigma, real alpha){
    real y = alpha/t;
    return log(alpha) - 2*log(t) + positive_normal_lpdf(y | nu, sigma);
  }
  
  // cdf of a t of a single accumulator
  real later_single_lcdf(real t, real nu, real sigma, real alpha){
    real y = alpha/t;
    return positive_normal_lccdf(y | nu, sigma);
  }
  
  // ccdf of a t of a single accumulator
  real later_single_lccdf(real t, real nu, real sigma, real alpha){
    real y = alpha/t;
    return positive_normal_lcdf(y | nu, sigma);
  }
  
  /*
  MULTIPLE ACCUMULATORS FUNCTIONS
  */
  // rng of rt and response given a bunch of accumulators
  vector later_rng(vector nu, vector sigma, vector alpha, real t0){
    int N_acc = num_elements(nu);
    vector[N_acc] rt = rep_vector(t0, N_acc);
    vector[2] output; // (rt, response)
    
    for(acc in 1:N_acc) rt[acc] += later_single_rng(nu[acc], sigma[acc], alpha[acc]);

    output[1] = min(rt);
    
    for(acc in 1:N_acc){
      if(rt[acc] == output[1]){
        output[2] = acc;
        break;
      }
    }
    
    return output;
  }
  
  // pdf of rt and response given a bunch of accumulators 
  real later_lpdf(real rt, int response, vector nu, vector sigma, vector alpha, real t0){
    real t = rt - t0;
    real log_lik = 0;
    int N_acc = num_elements(nu);
    
    for(acc in 1:N_acc){
      if(acc == response)
        log_lik += later_single_lpdf(t  | nu[acc], sigma[acc], alpha[acc]);
      else
        log_lik += later_single_lccdf(t | nu[acc], sigma[acc], alpha[acc]); 
    }
    
    return log_lik;
  }
  
  // function to convert real to integer
  int response_to_integer(real response, int N_acc){
    int output;

    for(i in 1:N_acc){
      if(response == i){
        output = i;
        break;
      }
    }

    return output;
  }
// }
