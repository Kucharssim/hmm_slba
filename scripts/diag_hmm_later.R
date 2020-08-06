args <- commandArgs(trailingOnly = TRUE)
subject <- as.character(args[1])

library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.24.0-rc1/")
cmdstan_version()
library(here)
library(bayesplot)
library(posterior)
library(tidyverse)
library(rstan)

# get functions to plot analytic LATER curves ----
positive_normal <- readLines(here("stan", "helpers", "positive_normal.stan"))
later_functions <- readLines(here("stan", "later", "later_functions.stan"))
all_functions <- sprintf("functions{\n%s\n%s\n}", paste(positive_normal, collapse = "\n"), paste(later_functions, collapse = "\n"))

expose_stan_functions(stanmodel = stan_model(model_code=all_functions))

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

curve(dlater(rt=x, 1, c(0.6, 0.4), c(0.2, 0.2), c(0.5, 0.5), 0.1), from = 0, to = 3)

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

curve(plater(rt=x, 1, c(0.7, 0.3), c(0.5, 0.5), c(0.5, 0.5), 0.1), from = 0, to = 3)

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

state_prob <- matrix(samples$summary(variables = "state_prob")$mean, ncol = 2, byrow = TRUE)


#### Plot results -----
cum_check <- dplyr::tibble(
  rt = data$rt,
  accumulator = data$accumulator,
  p11 = NA, p12 = NA, p21 = NA, p22 = NA
)

for(i in seq_len(nrow(cum_check))){
  if(cum_check$accumulator[i] == 1){
    cum_check$p11[i] <- state_prob[i,1]
    cum_check$p21[i] <- state_prob[i,2]
    cum_check$p12[i] <- 0
    cum_check$p22[i] <- 0
  } else{
    cum_check$p11[i] <- 0
    cum_check$p21[i] <- 0
    cum_check$p12[i] <- state_prob[i,1]
    cum_check$p22[i] <- state_prob[i,2]
  }
}

cum_check <- arrange(cum_check, rt)

cum_check <- mutate_at(cum_check, .vars = vars(starts_with("p")), .funs = list(cumsum))

# cum_check$p11 <- cum_check$p11 / sum(state_prob[,1])
# cum_check$p12 <- cum_check$p12 / sum(state_prob[,1])
# cum_check$p21 <- cum_check$p21 / sum(state_prob[,2])
# cum_check$p22 <- cum_check$p22 / sum(state_prob[,2])

p_cum_check <- cum_check %>%
  pivot_longer(cols = starts_with("p"), 
               names_to = c("state", "response"), names_pattern = "p(.)(.)",
               values_to = "probability") %>%
  mutate(state = sprintf("State %s", state), 
         response = ifelse(response == 1, "Correct", "Incorrect"),
         probability = probability/nrow(cum_check)) %>%
  ggplot(aes(x = rt, y = probability)) + #, group = response, col = response)) +
  geom_line(size = 1, colour = "red") +
  facet_wrap(state~response) +
  #ylim(0:1) + 
  ylab("Cumulative probability") + xlab("Time (sec)") +
  theme_bw(base_size = 20)

add_prediction <- function(iteration = 1, chain = 1){
  pars <- list()
  # pars$nu    <- matrix(samples$summary(variables = "nu_vec")$mean,    ncol = 2)
  # pars$sigma <- matrix(samples$summary(variables = "sigma_vec")$mean, ncol = 2)
  # pars$alpha <- matrix(samples$summary(variables = "alpha_vec")$mean, ncol = 2)
  # pars$t0    <- samples$summary(variables = "t0_vec")$mean
  
  pars$nu    <- matrix(as.vector(samples$draws(variables = "nu_vec")   [iteration,chain,]), ncol = 2)
  pars$sigma <- matrix(as.vector(samples$draws(variables = "sigma_vec")[iteration,chain,]), ncol = 2)
  pars$alpha <- matrix(as.vector(samples$draws(variables = "alpha_vec")[iteration,chain,]), ncol = 2)
  pars$t0    <- as.vector(samples$draws(variables = "t0_vec")[iteration,chain,])
  #browser()
  cum_check_2 <- expand_grid(rt = seq(min(cum_check$rt), max(cum_check$rt), length.out = 100),
                             response = 1:2,
                             state = 1:2,
                             probability = NA)
  for(i in seq_len(nrow(cum_check_2))) {
    cum_check_2$probability[i] <- plater(rt = cum_check_2$rt[i], response = cum_check_2$response[i], 
                                         nu = pars$nu[cum_check_2$state[i],], 
                                         sigma = pars$sigma[cum_check_2$state[i],],
                                         alpha = pars$alpha[cum_check_2$state[i],],
                                         t0    = pars$t0[cum_check_2$state[i]]) * mean( state_prob[,cum_check_2$state[i]] )
  }
  
  geom_line(size = 0.5, alpha = 0.1,
            data = cum_check_2 %>% mutate(state = sprintf("State %s", state), response = ifelse(response == 1, "Correct", "Incorrect")))
  
}

p_cum_check + 
  sapply(1:50, add_prediction) +
  geom_line(size = 2, colour = "red", alpha = 0.8)
