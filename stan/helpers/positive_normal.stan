// functions{
  // rng function
  real positive_normal_rng(real mu, real sigma){
    real p1 = normal_cdf(0, mu, sigma);  // cdf with lower bound
    real u  = uniform_rng(p1, 1);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf
  }
  
  // log density
  real positive_normal_lpdf(real y, real mu, real sigma){
    real lnominator   = normal_lpdf(y  | mu, sigma);
    real ldenominator = normal_lccdf(0 | mu, sigma);
    
    return lnominator - ldenominator;
  }
  
  // log cumulative probability
  real positive_normal_lcdf(real y, real mu, real sigma){
    real lnominator  = log_diff_exp(
      normal_lcdf(y  | mu, sigma),
      normal_lcdf(0 | mu, sigma)
    );
    
    real ldenominator = normal_lccdf(0 | mu, sigma);
    
    return lnominator - ldenominator;
  }
  
  // log complement of cumulative probability
  real positive_normal_lccdf(real y, real mu, real sigma){
    real lnominator   = normal_lccdf(y | mu, sigma);
    real ldenominator = normal_lccdf(0 | mu, sigma); 
    
    return lnominator - ldenominator;
  }
// }

