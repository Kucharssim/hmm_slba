library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.24.0-rc1/")
cmdstan_version()
library(here)
library(bayesplot)

hmm_later_prior_pred <- cmdstan_model(stan_file = here("stan", "later", "hmm_later_prior_pred.stan"), include_paths = here()) 

hyperparams <- list(
  # General info
  N_obs           = 200,
  N_acc           = 2,
  N_sts           = 2,
  
  # Design: 
  nu_indices      = matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE),
  sigma_indices   = matrix(1, nrow = 2, ncol = 2, byrow = TRUE),
  alpha_indices   = matrix(c(1,1, 2,2), nrow = 2, ncol = 2, byrow = TRUE),
  t0_indices      = rep(1, 2),
  
  # Hyperparameters:
  nu_alpha        = c(16, 4, 10, 10), # correct state + guessing state
  
  sigma_mu        = as.array(0.5),
  sigma_sigma     = as.array(0.05),
  
  alpha_mu        = c(1, 0.5),
  alpha_sigma     = c(0.2, 0.1),
  
  t0_beta         = as.array(5),
  
  init_prob_alpha = rep(5, 2),
  tran_prob_alpha = list(c(8, 2), c(2, 8))
)

generated_data <- hmm_later_prior_pred$sample(data = hyperparams, iter_warmup = 0, iter_sampling = 1000, fixed_param = TRUE)

generated_data$summary(variables = c("prop_state"))
generated_data$summary(variables = c("prop_responses"))
generated_data$summary(variables = c("mean_rt"))

mcmc_hist(generated_data$draws(variables = "prop_state"),     facet_args = list(scales = "fixed", dir = "v"))
mcmc_hist(generated_data$draws(variables = "prop_responses"), facet_args = list(scales = "fixed", dir = "v"))
mcmc_hist(generated_data$draws(variables = "mean_rt"),        facet_args = list(scales = "fixed", dir = "v"))


hist(as.vector(generated_data$draws(variables = "rt")[5,1,]), breaks=20)
# library(rstan)
# foo <- stan_model(file = here("stan", "later", "hmm_later_prior_pred.stan"))
# dat <- sampling(foo, hyperparams, algorithm = "Fixed_param", chains = 1)
