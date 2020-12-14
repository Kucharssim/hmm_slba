# prerequisites: 01_prior_predictives was run
# the simulation iteration nr. x (x between 1 and 1000) can be run from terminal as 
# Rscript --vanilla scripts/03_full_bayes_estimation.R x
# the bash script "run_fit_sbc.sh" runs all simulation iterations in succession
args <- commandArgs(trailingOnly = TRUE)
iter <- as.integer(args[1])


if(!is.na(iter)) cat("iteration:", iter, "\n") else stop("invalid iteration number\n")

library(cmdstanr)
set_cmdstan_path(readRDS("path_to_cmdstan.Rds"))
cmdstan_version()
library(here)

# stan models that are used in this script
hmm_later_prior_pred <- cmdstan_model(stan_file = here("stan", "later", "hmm_later_prior_pred.stan"), include_paths = here()) 
hmm_later <- cmdstan_model(stan_file = here("stan", "later", "hmm_later.stan"), include_paths = here())

# load the prior predictive data
generated_data <- readRDS(here("saves", "prior_predictives.Rds"))


hyperparams <- readRDS(here("saves", "hyperparams.Rds"))
stan_data   <- hyperparams
stan_data$rt <- as.vector(generated_data$draws("rt")[iter,])
stan_data$responses <- as.vector(generated_data$draws("responses")[iter,])

# function to generate starting values from the priors
initFun <- function(){
  prior <- hmm_later_prior_pred$sample(data = hyperparams, iter_warmup = 0, iter_sampling = 1, fixed_param = TRUE)
  out <- list(
    #init_prob = as.vector(prior$draws(variables = "init_prob")),
    
    #nu = round(as.vector(prior$draws(variables = "nu")), 2),
    sigma = as.vector(prior$draws(variables = "sigma")),
    alpha = as.vector(prior$draws(variables = "alpha")),
    t0 = as.vector(prior$draws(variables = "t0"))
  )
  
  out 
}


quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# mcmc fitting
# 1) fit the model 
# 2) check if model label switched (against true states), refit the model if it did
# 3) repeat at maximum five times, then terminate, otherwise save the results into a separate file 
mcmc <- function(iter, stan_data, true_states, max_tries = 5){
  finished <- FALSE
  i <- 1
  while(i <= max_tries && !finished){
    samples <- suppressMessages(suppressWarnings(quiet(
      hmm_later$sample(data = stan_data, chains = 1, parallel_chains = 1, 
                       iter_warmup = 500, iter_sampling = 1000, 
                       init = initFun, refresh = 0, show_messages = FALSE)
    )))
    finished <- length(samples$output_files(include_failed = FALSE)) != 0  
    if(finished){
      # assess label switching
      est_states <- apply(matrix(samples$summary("state_prob", mean)$mean, ncol = 2, byrow = TRUE), 1, which.max)
      p_agree <- mean(est_states == true_states)
      if(p_agree < 0.5) converged <- FALSE
    } else{
      rm(samples)
    }
    i <- i + 1
  }
  
  if(finished){
    #return(as_tibble(as_draws_df(samples$draws())))
    #return(samples)
    samples$save_object(file = here("saves", "sbc", sprintf("sim_%s.Rds", iter)))
  }
}

mcmc(iter, stan_data, as.vector(generated_data$draws("state")[iter,,]))
