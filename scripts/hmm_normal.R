library(here)
library(cmdstanr)
library(posterior)
library(tidyverse)
library(patchwork)
cmdstan_version()
set_cmdstan_path("/Users/skuchar/.cmdstan/cmdstan")
set_cmdstan_path("/Users/skuchar/.cmdstan/cmdstan-2.24.0-rc1/")
cmdstan_version()

model <- cmdstanr::cmdstan_model(here::here("stan", "normal.stan"))
hmm_normal_rng <- cmdstanr::cmdstan_model(here::here("stan", "hmm_normal_rng.stan"))
hmm_normal  <- cmdstanr::cmdstan_model(here::here("stan", "hmm_normal.stan"), include_paths=here::here("stan"))
hmm_normal2 <- cmdstanr::cmdstan_model(here::here("stan", "hmm_normal2.stan"))


specs <- list(N_obs = 1000, N_sts = 3)
specs$mu <- seq_len(specs$N_sts)*10
specs$sigma <- seq_len(specs$N_sts)
specs$state_init_prob <- rep(1/specs$N_sts, specs$N_sts)
specs$state_tran_prob <- 5*diag(specs$N_sts) + matrix(runif(specs$N_sts^2), ncol = specs$N_sts)
specs$state_tran_prob <- sweep(specs$state_tran_prob, 1, rowSums(specs$state_tran_prob), "/")
specs

data <- hmm_normal_rng$sample(data = specs, chains = 1, iter_warmup = 0, iter_sampling = 1, fixed_param = TRUE)

hmm_data <- list(
  N_obs = specs$N_obs,
  N_sts = specs$N_sts,
  y     = as.vector(data$draws(variables = "y"))
)
plot(hmm_data$y, type="b", pch = 19, col = as.vector(data$draws(variables = "state")))

smpl <- hmm_normal$sample(data = hmm_data, chains = 4, parallel_chains = 4, refresh = 250)
samples <- hmm_normal2$sample(data = hmm_data, chains = 4, parallel_chains = 4, refresh = 250)
samples$cmdstan_diagnose()
samples$summary(variables = c("mu", "sigma", "state_tran_prob"))

state_prob <- as_draws_df( samples$draws(variables = "state_prob") )
state_prob <- summarise_draws(state_prob, "mean", ~quantile2(.x))
state_prob <- state_prob %>% 
  dplyr::mutate(variable = gsub("state_prob", "", variable)) %>%
  dplyr::mutate(variable = gsub("\\[", "", variable)) %>%
  dplyr::mutate(variable = gsub("\\]", "", variable)) %>%
  tidyr::separate(col="variable", into=c("state", "time"), sep=",", convert=TRUE) %>%
  dplyr::mutate(state = as.factor(state))
  
dplyr::tibble(
  time  = seq_len(hmm_data$N_obs),
  state = as.factor(as.vector(data$draws(variables = "state"))),
  y     = as.vector(data$draws(variables = "y"))
) %>%
  ggplot2::ggplot(ggplot2::aes(x=time, y=y)) +
  ggplot2::geom_line() + 
  ggplot2::geom_point(ggplot2::aes(col=state), size = 2) + 
state_prob %>% 
  ggplot2::ggplot(ggplot2::aes(x=time, y=mean, group=state, fill=state, col=state)) +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(ggplot2::aes(ymin=q5, ymax=q95), alpha = 0.5, lwd=0) +
  patchwork::plot_layout(nrow = 2)


library(rstan)
rstan_model <- rstan::stan_model(here::here("stan", "hmm_normal.stan"))
rstan_fit   <- rstan::sampling(rstan_model, hmm_data, cores=4,chains=4)

stan_data <- list(
  N_obs = 1000,
  N_sts = 3,
  mu    = c(0, 5, 15),
  sigma = c(3, 1, 7),
  
  state_init_prob = c(0.3, 0.6, 0.1),
  state_tran_prob = list(c(0.8, 0.05, 0.15), c(0.3, 0.4, 0.3), c(0.1, 0.1, 0.8))
)

generated_data <- hmm_normal_rng$sample(data = stan_data, chains = 1, iter_warmup = 0, iter_sampling = 1, fixed_param = TRUE)

stan_data  <- list(
  N_obs = stan_data$N_obs,
  N_sts = stan_data$N_sts,
  y     = as.vector(generated_data$draws(variables = "y"))
)

stan_fit <- sampling(rstan_model, stan_data, cores = 4, chains = 4)
cmdstan_fit2 <- hmm_normal2$sample(stan_data, chains = 4, parallel_chains = 4, refresh = 250)
cmdstan_fit2$summary(variables = c("mu", "sigma", "state_tran_prob"))
