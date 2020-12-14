library(cmdstanr)
set_cmdstan_path(readRDS("path_to_cmdstan.Rds"))
cmdstan_version()
library(here)
library(bayesplot)
library(posterior)
library(tidyverse)

### Prior predictives ----
n_predictive <- 1000
# stan model for generating data from the model, given priors on the parameters
hmm_later_prior_pred <- cmdstan_model(stan_file = here("stan", "later", "hmm_later_prior_pred.stan"), include_paths = here()) 

hyperparams <- list(
  # General info
  N_obs           = 200, # number of observations per sim
  N_acc           = 2,   # number of accumulators (response options)
  N_sts           = 2,   # number of latent states
  
  # Design: 4 drift rates, 1 std of drift rates, 2 boundaries (each for each state), one non-decision time
  nu_indices      = matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE),
  sigma_indices   = matrix(1, nrow = 2, ncol = 2, byrow = TRUE),
  alpha_indices   = matrix(c(1,1, 2,2), nrow = 2, ncol = 2, byrow = TRUE),
  t0_indices      = rep(1, 2),
  
  # Hyperparameters:
  nu_alpha        = c(14, 6, 10, 10), # correct state + guessing state
  
  sigma_mu        = as.array(0.4),
  sigma_sigma     = as.array(0.1),
  
  alpha_mu        = c(0.5, 0.25),
  alpha_sigma     = c(0.1, 0.05),
  
  t0_beta         = as.array(5),
  
  init_prob_alpha = rep(5, 2),
  tran_prob_alpha = list(c(8, 2), c(2, 8))
)
saveRDS(hyperparams, here("saves", "hyperparams.Rds")) # save the hyperparameters (so that we can reuse for fitting)

# generate data (change seed to change results)
generated_data <- hmm_later_prior_pred$sample(data = hyperparams, iter_warmup = 0, iter_sampling = n_predictive, fixed_param = TRUE, seed = 2020)
generated_data$save_object(here("saves", "prior_predictives.Rds"))
generated_data <- readRDS(here("saves", "prior_predictives.Rds"))

# extract variables into separate objects for better manipulation
rt <- generated_data$draws("rt")
responses <- generated_data$draws("responses")
state <- generated_data$draws("state")

# average run length
run_lengths <- sapply(1:dim(state)[1], function(x) mean(rle(as.vector(state[x,,]))$lengths))
summary(run_lengths)
sd(run_lengths)
#

rt_df <- as_tibble(as_draws_df(rt)) %>%
  pivot_longer(cols = starts_with("rt"), names_to = "trial", values_to = "rt", names_pattern = "rt\\[(.*)\\]")

responses_df <- as_tibble(as_draws_df(responses)) %>%
  pivot_longer(cols = starts_with("responses"), names_to = "trial", values_to = "response", names_pattern = "responses\\[(.*)\\]")

state_df <- as_tibble(as_draws_df(state)) %>%
  pivot_longer(cols = starts_with("state"), names_to = "trial", values_to = "state", names_pattern = "state\\[(.*)\\]")

# generate descriptives of the prior predictive response times, under the two states, and for the two alternatives
prior_pred_df <- left_join(rt_df, responses_df) %>% left_join(state_df) %>%
  mutate(state_text = ifelse(state == 1, "Controlled", "Guessing"),
         response_text = ifelse(response == 1, "Correct", "Error")) %>%
  mutate(comb_state_response = sprintf("State: %s, Response: %s", state_text, response_text))

prior_pred_df %>% group_by(state_text, response_text) %>% summarise(long_rt = mean(rt > 3))
prior_pred_df %>% group_by(.iteration) %>% summarise(long_rt = mean(rt > 3)) %>% ungroup() %>% summary()

# descriptives: RT per state and response
prior_pred_df %>%
  group_by(.draw, state_text, response_text) %>%
  summarise(mean_rt = mean(rt)) %>%
  ungroup() %>%
  group_by(state_text, response_text) %>%
  summarise(mean  = mean(mean_rt), 
            sd    = sd(mean_rt),
            lower = quantile(mean_rt, 0.025),
            lqrt  = quantile(mean_rt, 0.25),
            median= quantile(mean_rt, 0.5),
            uqrt  = quantile(mean_rt, 0.75),
            upper = quantile(mean_rt, 0.975)) %>%
  knitr::kable(digits = 2)
# knitr::kable(digits = 2, format = "latex", booktabs = TRUE)

prior_pred_df %>%
  group_by(.draw, comb_state_response) %>%
  summarise(rt = mean(rt)) %>%
  ungroup() %>%
  ggplot(aes(x = rt)) + 
  geom_histogram(bins = 50) +
  facet_wrap(~comb_state_response, scales = "free_y") + 
  ylab("Count") + xlab("Average response time (sec)") +
  theme_bw(base_size = 15) + 
  theme(aspect.ratio = 0.7, 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())
ggsave(filename = here("figures", "simulations", "prior_pred_mean_rt.png"), width = 20, height = 15, unit = "cm")

prior_pred_df %>% 
  subset(.draw < 10 & rt < 3) %>%
  ggplot(aes(x = rt, group = .draw)) + 
  geom_density() + 
  facet_wrap(~comb_state_response, scales = "free_y") + 
  ylab("Density") + xlab("Response time (sec)") +
  theme_bw(base_size = 15) + 
  theme(aspect.ratio = 0.7, 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())
ggsave(filename = here("figures", "simulations", "prior_pred_rt.png"), width = 20, height = 15, unit = "cm")

# descriptives: responses per state
prior_pred_df %>%
  group_by(.draw, state_text) %>%
  summarise(prop_correct = mean(response == 1), n_correct = sum(response == 1)) %>%
  ungroup() %>%
  group_by(state_text) %>% 
  summarise(mean  = mean(prop_correct), 
            sd    = sd(prop_correct),
            lower = quantile(prop_correct, 0.025),
            lqrt  = quantile(prop_correct, 0.25),
            median= quantile(prop_correct, 0.5),
            uqrt  = quantile(prop_correct, 0.75),
            upper = quantile(prop_correct, 0.975)) %>%
  knitr::kable(digits = 2)
# knitr::kable(digits = 2, format = "latex", booktabs = TRUE)

prior_pred_df %>% 
  group_by(.draw, state_text) %>%
  summarise(prop_correct = mean(response == 1)) %>%
  ggplot(aes(x = prop_correct)) +
  geom_histogram(bins = 50) +
  facet_wrap(~state_text) + 
  ylab("Count") + xlab("Proportion of correct responses") +
  theme_bw(base_size = 15) + 
  theme(aspect.ratio = 0.7, 
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())
ggsave(filename = here("figures", "simulations", "prior_pred_accuracy.png"), width = 20, height = 8, unit = "cm")
