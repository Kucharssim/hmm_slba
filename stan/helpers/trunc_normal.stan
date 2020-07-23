// functions{
  // rng function
  real trunc_normal_rng(real mu, real sigma, real lb, real ub){
    real p1 = normal_cdf(lb, mu, sigma);  // cdf with lower bound
    real p2 = normal_cdf(ub, mu, sigma);  // cdf with upper bound
    real u = uniform_rng(p1, p2);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf
  }
  
  // log density
  real trunc_normal_lpdf(real y, real mu, real sigma, real lb, real ub){
    real lnominator   = normal_lpdf(y | mu, sigma);
    real ldenominator = log_diff_exp(
      normal_lcdf(ub | mu, sigma),
      normal_lcdf(lb | mu, sigma)
    );
    
    return lnominator - ldenominator;
  }
  
  // log cumulative probability
  real trunc_normal_lcdf(real y, real mu, real sigma, real lb, real ub){
    real lnominator  = log_diff_exp(
      normal_lcdf(y  | mu, sigma),
      normal_lcdf(lb | mu, sigma)
    );
    
    real ldenominator = log_diff_exp(
      normal_lcdf(ub | mu, sigma),
      normal_lcdf(lb | mu, sigma)
    );
    
    return lnominator - ldenominator;
  }
  
  // log complement of cumulative probability
  real trunc_normal_lccdf(real y, real mu, real sigma, real lb, real ub){
    real lnominator  = log_diff_exp(
      normal_lcdf(ub | mu, sigma),
      normal_lcdf(y  | mu, sigma)
    );
    
    real ldenominator = log_diff_exp(
      normal_lcdf(ub | mu, sigma),
      normal_lcdf(lb | mu, sigma)
    );
    
    return lnominator - ldenominator;
  }
// }

