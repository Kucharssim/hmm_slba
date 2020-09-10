library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.24.0-rc1/")
cmdstan_version()
library(here)
library(bayesplot)
library(posterior)
library(tidyverse)
library(patchwork)

source(here("R", "expose_stan_functions.R"))

# compute the cumulative rt distributions  and specific draw from the parameters (samples_params)
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
  # if(accumulator == 2) n_cases <- nrow(data) - n_cases
  plater(rt, accumulator, nu, sigma, alpha, t0) * n_cases / nrow(data)
}
# compute the pdf for each state and specific draw from the parameters (samples_params)
call_dlater <- function(rt, accumulator, draw_, state_){
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
  #if(accumulator == 2) n_cases <- nrow(data) - n_cases
  dlater(rt, accumulator, nu, sigma, alpha, t0) * n_cases / nrow(data)
}

# posterior prediction of the later model for each observation and a specific draw from the parameters (samples_params)
call_rlater <- function(state_, draw_){
  #obs_state_prob <- subset(state_prob, .draw == draw_ & obs == obs_)$probability 
  #state_ <- sample(1:2, 1, TRUE, obs_state_prob)
  nu <- subset(samples_params, .draw == draw_ & parameter == "nu_acc" & state == state_)$value
  nu <- c(nu, 1-nu)
  sigma <- subset(samples_params, .draw == draw_ & parameter == "sigma")$value
  sigma <- rep(sigma, 2)
  alpha <- subset(samples_params, .draw == draw_ & parameter == "alpha" & state == state_)$value
  alpha <- rep(alpha, 2)
  t0    <- subset(samples_params, .draw == draw_ & parameter == "t0")$value

  rlater(1, nu, sigma, alpha, t0)
}

#### Read data and models -----
hmm_later_prior_pred <- cmdstan_model(stan_file = here("stan", "later", "hmm_later_prior_pred.stan"), include_paths = here()) 
hmm_later <- cmdstan_model(stan_file = here("stan", "later", "hmm_later.stan"), include_paths = here())
parameters       <- c("nu_vec[1,1]", "nu_vec[2,1]",  "alpha",  "sigma",   "t0", "init_prob[1]", "tran_prob[1,1]", "tran_prob[2,2]")
parameters_latex <- c("\\nu_1^{(1)}", "\\nu_1^{(2)}", "\\alpha^{(1)}", "\\alpha^{(2)}", "\\sigma", "\\tau", "\\pi_1",        "\\rho_{11}",      "\\rho_{22}")
parameters_latex <- sprintf("$%s$", parameters_latex)
samples_list <- list()

for(subject in LETTERS[1:11]){
  cat("===== Subject", subject, "=====\n")
  data <- readr::read_csv(here::here("data", sprintf("dutilh_2010_subject_%s.csv", subject)), 
                          col_types = cols(
                            pacc = col_double(),
                            rt = col_double(),
                            accumulator = col_integer()
                          ))
  data$obs <- seq_len(nrow(data))
  
  #### Get MCMC samples
  samples <- readRDS(here("saves", "fit_hmm_later", sprintf("dutilh_2010_subject_%s.Rds", subject)))
  samples_list[[subject]] <- samples
  samples$summary(variables = parameters) %>% print()
  
  traceplot <- bayesplot::mcmc_trace(samples$draws(variables = parameters), facet_args = list(ncol = 3)) +
    ggplot2::ggtitle(sprintf("Participant %s", subject))
  ggplot2::ggsave(filename = here("figures", "traceplots", sprintf("dutilh_2010_subject_%s.png", subject)), 
                  width = 15, height = 15, units = "cm")
  
  # reshape state probability draws into a long data frame
  state_prob <- samples$draws("state_prob") %>% 
    as_draws_df() %>% 
    pivot_longer(cols          = starts_with("state_prob"), 
                 names_pattern = "state_prob\\[(.*),(.*)\\]", 
                 names_to      = c("state", "obs"),
                 values_to     = "probability")  %>%
    mutate(state = as.integer(state), obs = as.integer(obs))
  
  # descriptive plot
  state_prob_summary <- state_prob %>% 
    subset(state == 1) %>% 
    group_by(obs) %>% 
    summarise(mean = mean(probability), lower = quantile(probability, 0.05), upper = quantile(probability, 0.95))
  desc_plot <- list()
  desc_plot[['response']] <- ggplot(data, aes(x = obs, y = accumulator)) + 
    geom_line(data = data %>% mutate(ma = forecast::ma(accumulator, 5)), mapping = aes(x = obs, y = ma)) + 
    geom_point() + 
    xlab("") + ylab("Response") +
    scale_y_continuous(breaks = 1:2, labels = c("Incorrect", "Correct"), expand = c(0.2,0)) +
    theme_light()
  desc_plot[['rt']] <- ggplot(data, aes(x = obs, y = rt)) + 
    geom_line() +
    xlab("") + ylab("Response time (sec)") +
    scale_y_continuous(limits = c(0, NA)) +
    theme_light()
  desc_plot[['state']] <- ggplot(state_prob_summary, aes(x = obs, y = mean, ymin = lower, ymax = upper)) + 
    geom_line() + 
    geom_ribbon() +
    xlab("Trial") + ylab("P(controlled state)") + 
    theme_light()
  desc_plot[['response']] + desc_plot[['rt']] + desc_plot[['state']] + patchwork::plot_layout(nrow = 3) +
    patchwork::plot_annotation(title = sprintf("Participant %s", subject)) &
    theme(text = element_text(size = 18), legend.title = element_blank(), legend.position = "none")
  ggplot2::ggsave(filename = here("figures", "desc_plots", sprintf("dutilh_2010_subject_%s.png", subject)), 
                  width = 10, height = 20, units = "cm")
  
  random_draws <- sample(unique(state_prob$.draw), size = 1000, replace = FALSE)
  
  # get samples of parameters
  samples_params <- samples$draws(parameters) %>%
    as_draws_df() %>%
    as_tibble() %>%
    pivot_longer(cols = starts_with(parameters), names_pattern = "(.*)\\[(.*)\\]", names_to = c("parameter", "state")) %>%
    subset(.draw %in% random_draws) %>%
    mutate(parameter = gsub("_vec", "_acc", parameter), state = as.integer(gsub(",(1|2)", "", state)))
  
  
  # posterior predictives
  #post_pred <- expand.grid(.draw = random_draws[1:100], .obs = seq_len(nrow(data)), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  post_pred <- expand.grid(.draw = random_draws, state = 1:2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  post_pred$rt <- numeric(length = nrow(post_pred))
  post_pred$accumulator <- integer(length = nrow(post_pred))
  for(i in seq_len(nrow(post_pred))){
    pp                       <- call_rlater(post_pred$state[i], post_pred$.draw[i])
    post_pred$rt[i]          <- pp[[1]]
    post_pred$accumulator[i] <- pp[[2]]
    if(i %% 100 == 0) cat("iter", i, "done\n")
  }
  
  post_pred <- state_prob %>% 
    subset(.draw %in% random_draws) %>% 
    left_join(post_pred) %>%
    group_by(.draw, obs) %>% 
    mutate(which_state = sample(1:2, 1, TRUE, probability)) %>%
    ungroup() %>%
    subset(state == which_state)
  
  post_pred_rt_plot <- post_pred %>%
    group_by(obs) %>%
    summarise(median = median(rt), 
              lower80 = quantile(rt, 0.1),  upper80 = quantile(rt, 0.9),
              lower50 = quantile(rt, 0.25), upper50 = quantile(rt, 0.75)) %>%
    ggplot(aes(x=obs, y=median)) +
      geom_ribbon(mapping = aes(x=obs, ymin=lower80, ymax=upper80), inherit.aes = FALSE, fill = "#DCBCBC", alpha = 0.7) +
      geom_ribbon(mapping = aes(x=obs, ymin=lower50, ymax=upper50), inherit.aes = FALSE, fill = "#B97C7C", alpha = 0.7) +
      geom_line(col = "#8F2727", size = 1) + 
      geom_line(aes(x=obs, y=rt), data = data, inherit.aes = FALSE, col = "black", size = 0.75) +
      xlab("Trial") + ylab("Response time (sec)") + ggtitle(sprintf("Participant %s", subject)) +
      theme_light() +
      theme(text = element_text(size = 15), legend.title = element_blank(), legend.position = "none")
  ggplot2::ggsave(plot = post_pred_rt_plot, 
                  filename = here("figures", "post_pred_rt", sprintf("dutilh_2010_subject_%s.png", subject)),
                  width = 20, height = 7.5, units = "cm")
  
  post_pred_response_plot <- post_pred %>%
    group_by(obs) %>%
    summarise(mean = mean(accumulator == 1)) %>%
    ungroup() %>%
    ggplot(aes(x=obs,y=mean)) +
    geom_line(col = "#8F2727", size = 1) +
    geom_line(mapping = aes(x = obs, y = ma), data = data %>% mutate(ma = forecast::ma(ifelse(accumulator==1, 1, 0), 10)), inherit.aes = FALSE, size = 1) + 
    geom_point(mapping = aes(x=obs, y=ifelse(accumulator==1, 1, 0)), data = data, inherit.aes = FALSE, size = 0.5) +
    xlab("Trial") + ylab("Response") + ggtitle(sprintf("Participant %s", subject)) +
    ylim(0, 1) +
    theme_light() +
    theme(text = element_text(size = 15), legend.title = element_blank(), legend.position = "none")
  ggplot2::ggsave(plot = post_pred_response_plot, 
                  filename = here("figures", "post_pred_response", sprintf("dutilh_2010_subject_%s.png", subject)),
                  width = 20, height = 7.5, units = "cm")
  
  # plot empirical cumulative rt distributions scaled by proportions of correct/incorrect and state1/state2
  cdf_plot <- state_prob %>%
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
    xlab("Response time (sec)") +
    ggtitle(sprintf("Participant %s", subject)) +
    theme_light() +
    theme(text = element_text(size = 20), legend.title = element_blank(), legend.position = "none")
  
  
  # # plot model based cumulative rt distributions
  df <- expand.grid(
    rt          = seq(0, max(data$rt), length.out = 20),
    accumulator = 1:2,
    state       = 1:2,
    .draw       = random_draws[1:100],
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)


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

  cdf_plot <- cdf_plot +
    geom_line(data = df, linetype = 2) +
    geom_ribbon(aes(ymin=lower,ymax=upper), alpha = 0.25, data = df)
  ggplot2::ggsave(plot = cdf_plot, filename = here("figures", "cdf_plots", sprintf("dutilh_2010_subject_%s.png", subject)),
                  width = 15, height = 15, units = "cm")
  
  pdf_plot <- data %>% 
    mutate(response = ifelse(accumulator == 1, "Correct", "Incorrect")) %>%
    ggplot(aes(x = rt)) + 
      geom_histogram(bins = 25, fill = "grey", col = "black") + 
      facet_wrap(~response)
  
  df <- expand.grid(
    rt          = seq(0, max(data$rt), length.out = 100),
    accumulator = 1:2,
    state       = 1:2,
    .draw       = random_draws[1:100],
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  
  
  df <- df %>%
    group_by(accumulator, state, .draw) %>%
    # compute the pdf
    mutate(pdf = call_dlater(rt, unique(accumulator), unique(.draw), unique(state))) %>%
    ungroup() %>%
    mutate(accumulator = ifelse(accumulator == 1, "Correct", "Incorrect"))
  df_by_state <- df %>%  
    # compute confidence bounds + median
    group_by(state, accumulator, rt) %>%
    summarise(median = median(pdf), lower = quantile(pdf, 0.05), upper = quantile(pdf, 0.95)) %>%
    ungroup()
  df_aggregate <- df %>% 
    pivot_wider(id_cols = c("rt", "accumulator", ".draw"), names_from = state, values_from = pdf) %>%
    mutate(pdf = `1` + `2`) %>%
    group_by(accumulator, rt) %>%
    summarise(median = median(pdf), lower = quantile(pdf, 0.05), upper = quantile(pdf, 0.95))
  
  pdf_plot <- data %>% 
    mutate(prop = ifelse(accumulator == 1, mean(data$accumulator == 1), mean(data$accumulator == 2))) %>%
    mutate(cut = cut(rt, breaks = seq(0, max(rt), length.out = 50))) %>%
    group_by(cut, accumulator) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    mutate(prop = ifelse(accumulator == 1, n / sum(accumulator == 1), n / sum(accumulator == 2))) %>%
    mutate(cut = gsub("\\(|\\]", "", cut)) %>%
    separate(cut, c("min", "max"), sep = ",", convert = TRUE) %>%
    mutate(x = (min + max) /2, width = max - min) %>%
    mutate(prop = (prop / sum(prop)) / width) %>%
    mutate(accumulator = ifelse(accumulator == 1, "Correct", "Incorrect")) %>%
    ggplot(aes(x=x, y=prop)) + geom_bar(stat="identity", fill = "grey", col = "black", alpha = 0.3) + 
      geom_ribbon(data = df_by_state, mapping  = aes(x = rt, ymin = lower, ymax = upper, group = state, fill = as.factor(state)), alpha = 0.5, inherit.aes = FALSE) + 
      geom_line(data   = df_by_state, mapping  = aes(x = rt, y = median, group = state), alpha = 0.5, inherit.aes = FALSE) + 
      geom_ribbon(data = df_aggregate, mapping = aes(x = rt, ymin = lower, ymax = upper), fill = "black", alpha = 0.7, inherit.aes = FALSE) +
      geom_line(data   = df_aggregate, mapping = aes(x = rt, y = median), inherit.aes = FALSE) +
      xlab("Response time (sec)") + ylab("Density") +
      facet_wrap(~accumulator) + 
      ggtitle(sprintf("Participant %s", subject)) +
      theme_light() +
      theme(text = element_text(size = 20), legend.title = element_blank(), legend.position = "none")
  ggplot2::ggsave(plot = pdf_plot, filename = here("figures", "pdf_plots", sprintf("dutilh_2010_subject_%s.png", subject)),
                  width = 25, height = 15, units = "cm")
  
}


samples_list[["K"]]$summary(variables = parameters) %>% 
  rename("Parameter"  = "variable",
         "Mean"       = "mean",
         "Median"     = "median",
         "SD"         = "sd",
         "5\\%"       = "q5",
         "95\\%"      = "q95",
         "$\\hat{R}$" = "rhat",
         "Bulk"       = "ess_bulk",
         "Tail"       = "ess_tail") %>%
  select(-mad) %>%
  mutate(Parameter=parameters_latex) %>%
  knitr::kable(digits  = c(0, 2, 2, 2, 2, 2, 3, 0, 0),
               format  = "latex", booktabs = TRUE, escape = FALSE, linesep = c('', '', '', '', '', '\\addlinespace'),
               caption = "Descriptives of the posterior draws for Participant K from \\citet{dutilh2011phase}.",
               label   = "tab:pars_K") %>%
  kableExtra::add_header_above(c(" " = 4, "Quantile" = 2, " " = 1, "ESS" = 2))
