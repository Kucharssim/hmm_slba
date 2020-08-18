library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.24.0-rc1/")
cmdstan_version()
library(here)
library(bayesplot)
library(posterior)
library(tidyverse)

source(here("R", "expose_stan_functions.R"))


#### Read data and models -----
hmm_later_prior_pred <- cmdstan_model(stan_file = here("stan", "later", "hmm_later_prior_pred.stan"), include_paths = here()) 
hmm_later <- cmdstan_model(stan_file = here("stan", "later", "hmm_later.stan"), include_paths = here())
parameters <- c("nu_vec[1,1]", "nu_vec[2,1]", "alpha", "sigma", "t0")

for(subject in LETTERS[1:11]){
  cat("===== Subject", subject, " =====\n")
  data <- readr::read_csv(here::here("data", sprintf("dutilh_2010_subject_%s.csv", subject)), 
                          col_types = cols(
                            pacc = col_double(),
                            rt = col_double(),
                            accumulator = col_integer()
                          ))
  data$obs <- seq_len(nrow(data))
  
  #### Get MCMC samples
  samples <- readRDS(here("saves", "fit_hmm_later", sprintf("dutilh_2010_subject_%s.Rds", subject)))
  samples$summary(variables = c("nu_vec[1,1]", "nu_vec[2,1]", "alpha", "sigma", "t0")) %>% print()
  
  traceplot <- bayesplot::mcmc_trace(samples$draws(variables = parameters), facet_args = list(ncol = 2)) +
    ggplot2::ggtitle(sprintf("Subject %s", subject))
  ggplot2::ggsave(filename = here("figures", "traceplots", sprintf("dutilh_2010_subject_%s.png", subject)))
  
  # reshape state probability draws into a long data frame
  state_prob <- samples$draws("state_prob") %>% 
    as_draws_df() %>% 
    pivot_longer(cols          = starts_with("state_prob"), 
                 names_pattern = "state_prob\\[(.*),(.*)\\]", 
                 names_to      = c("state", "obs"),
                 values_to     = "probability")  %>%
    mutate(state = as.integer(state), obs = as.integer(obs))
  
  random_draws <- sample(unique(state_prob$.draw), size = 1000, replace = FALSE)
  
  # plot empirical cumulative rt distributions scaled by proportions of correct/incorrect and state1/state2
  cum_plot <- state_prob %>%  
    # subset only a portion of draws
    subset(.draw %in% random_draws) %>%
    left_join(data) %>%
    arrange(.draw, state, accumulator, rt) %>%
    group_by(.draw, state, accumulator) %>%
    mutate(cum_prob = cumsum(probability)/nrow(data)) %>%
    ungroup() %>%
    group_by(state, accumulator, rt) %>%
    summarise(median = median(cum_prob), lower = quantile(cum_prob, 0.05), upper = quantile(cum_prob, 0.95)) %>%
    mutate(group = sprintf("%s state, %s response", 
                           ifelse(state==1, "accurate", "fast"),
                           ifelse(accumulator==1, "correct", "incorrect")
                           )
           ) %>%
    ggplot(aes(x = rt, y = median, group = group, fill = group)) + 
    geom_line() + 
    geom_ribbon(aes(ymin=lower,ymax=upper), alpha = 0.5) +
    scale_fill_manual(values = c("darkred", "red", "darkblue", "blue")) + 
    ylab("Cumulative probability") +
    xlab("Response time") + 
    ggtitle(sprintf("Subject %s", subject)) + 
    theme_light() + 
    theme(text = element_text(size = 20), legend.title = element_blank(), legend.position = "none")
  
  # get samples of parameters
  samples_params <- samples$draws(parameters) %>%
    as_draws_df() %>%
    as_tibble() %>%
    pivot_longer(cols = starts_with(parameters), names_pattern = "(.*)\\[(.*)\\]", names_to = c("parameter", "state")) %>%
    subset(.draw %in% random_draws) %>%
    mutate(parameter = gsub("_vec", "_acc", parameter), state = as.integer(gsub(",1", "", state)))
  
  
  # plot model based cumulative rt distributions
  df <- expand.grid(
    rt          = seq(0, max(data$rt), length.out = 20),
    accumulator = 1:2,
    state       = 1:2,
    .draw       = random_draws[1:400],
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  # compute the cumulative rt distributions
  call_plater <- function(rt, accumulator, draw_, state_){
    #browser()
    nu <- subset(samples_params, .draw == draw_ & parameter == "nu_acc" & state == state_)$value
    nu <- c(nu, 1-nu)
    sigma <- subset(samples_params, .draw == draw_ & parameter == "sigma")$value
    sigma <- rep(sigma, 2)
    alpha <- subset(samples_params, .draw == draw_ & parameter == "alpha" & state == state_)$value
    alpha <- rep(alpha, 2)
    t0    <- subset(samples_params, .draw == draw_ & parameter == "t0")$value
    
    # scale by the probability of response and state
    n_cases <- sum(subset(state_prob, .draw == draw_ & state == state_)$probability)
    if(accumulator == 2) n_cases <- nrow(data) - n_cases
    plater(rt, accumulator, nu, sigma, alpha, t0) * n_cases / nrow(data)
  }
  
  df <- df %>% 
    group_by(accumulator, state, .draw) %>%
    # compute the cdf
    mutate(cum_prob = call_plater(rt, unique(accumulator), unique(.draw), unique(state))) %>%
    ungroup() %>%
    # compute confidence bounds + median
    group_by(state, accumulator, rt) %>%
    summarise(median = median(cum_prob), lower = quantile(cum_prob, 0.05), upper = quantile(cum_prob, 0.95)) %>%
    ungroup() %>%
    # labeling should be the same as for the data plot
    mutate(group = sprintf("%s state, %s response", 
                           ifelse(state==1, "accurate", "fast"),
                           ifelse(accumulator==1, "correct", "incorrect")
                           )
           )
    
  cum_plot +
    geom_line(data = df, linetype = 2) + 
    geom_ribbon(aes(ymin=lower,ymax=upper), alpha = 0.25, data = df)
  ggplot2::ggsave(filename = here("figures", "cum_plots", sprintf("dutilh_2010_subject_%s.png", subject)))

}
