args <- commandArgs(trailingOnly = TRUE)
subject <- as.character(args[1])

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.24.0-rc1/")
cmdstan_version()
library(here)

#### Read data and models -----
data <- readr::read_csv(here::here("data", sprintf("dutilh_2010_subject_%s.csv", subject)))

hmm_later_prior_pred <- cmdstan_model(stan_file = here("stan", "later", "hmm_later_prior_pred.stan"), include_paths = here()) 
hmm_later <- cmdstan_model(stan_file = here("stan", "later", "hmm_later.stan"), include_paths = here())
hyperparams <- readRDS(here("saves", "hyperparams.Rds"))

stan_data <- hyperparams
stan_data$N_obs <- nrow(data)
stan_data$rt <- data$rt
stan_data$responses <- data$accumulator

initFun <- function(){
  prior <- hmm_later_prior_pred$sample(data = hyperparams, iter_warmup = 0, iter_sampling = 1, fixed_param = TRUE)
  out <- list(
    sigma = as.vector(prior$draws(variables = "sigma")),
    alpha = as.vector(prior$draws(variables = "alpha")),
    t0 = as.vector(prior$draws(variables = "t0"))
  )
  
  out 
}

#### Get MCMC samples ----
samples <- hmm_later$sample(data = stan_data, 
                            chains = 8, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                            init = initFun, adapt_delta = 0.99)
samples$summary(variables = c("nu_vec[1,1]", "nu_vec[2,1]", "alpha", "sigma", "t0"))

saveRDS(samples, file = here("saves", "fit_hmm_later", sprintf("dutilh_2010_subject_%s.Rds", subject)))
