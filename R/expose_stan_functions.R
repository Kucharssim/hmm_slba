# expose LATER functions from Stan
library(rstan)
library(here)
# reads functions definint truncated normal distribution (0, Inf)
positive_normal <- readLines(here("stan", "helpers", "positive_normal.stan"))
# reads functions defining the LATER distributions
later_functions <- readLines(here("stan", "later", "later_functions.stan"))

all_functions <- sprintf("functions{\n%s\n%s\n}", paste(positive_normal, collapse = "\n"), paste(later_functions, collapse = "\n"))

expose_stan_functions(stanmodel = stan_model(model_code=all_functions))
rm(all_functions, positive_normal, later_functions)


# make vectorized analogues to the functions which are similar to the R style of defining distributions
# pdf ---
dlater_single <- function(rt, nu, sigma, alpha, t0, log = FALSE){
  out <- sapply(rt, function(rt_){
    if (rt_ > t0) {
      later_single_lpdf(rt_ - t0, nu, sigma, alpha)
    } else {
      -Inf
    }
  })
  
  if (log) {
    return(out)
  } else {
    return(exp(out))
  }
}

dlater <- function(rt, response, nu, sigma, alpha, t0, log = FALSE){
  out <- sapply(rt, function(rt_){
    if(rt_ > t0) {
      later_lpdf(rt_, response, nu, sigma, alpha, t0)
    } else {
      -Inf
    }
  })
  
  if(log) {
    return(out)
  } else {
    return(exp(out))
  }
}

# cdf ---
plater_single <- function(rt, nu, sigma, alpha, t0, log.p = FALSE){
  out <- sapply(rt, function(rt_) {
    # log(integrate(dlater_single, lower = 0, upper = rt_, nu=nu, sigma=sigma, alpha=alpha, t0=t0)$value)
    if (rt_ > t0) {
      later_single_lcdf(rt_ - t0, nu, sigma, alpha)
    } else {
      -Inf
    }
  })
  
  if(log.p){
    return(out)
  } else{
    return(exp(out))
  }
}

plater <- function(rt, response, nu, sigma, alpha, t0, log.p = FALSE){
  out <- sapply(rt, function(rt_) {
    log(integrate(dlater, lower = 0, upper = rt_, response=response, nu=nu, sigma=sigma, alpha=alpha, t0=t0)$value)
  })
  
  if(log.p){
    return(out)
  } else{
    return(exp(out))
  }
}

# rng
rlater_single <- function(n, nu, sigma, alpha, t0) {
  replicate(n, later_single_rng(nu, sigma, alpha)) + t0
}

rlater <- function(n, nu, sigma, alpha, t0){
  out <- replicate(n, later_rng(nu, sigma, alpha, t0))
  out <- as.data.frame(t(out))
  colnames(out) <- c("rt", "response")
  
  out
}
