library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.24.0-rc1/")
cmdstan_version()
library(here)
library(bayesplot)
library(posterior)
library(tidyverse)

### Prior predictives ----
n_predictive <- 1000
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
  nu_alpha        = c(14, 6, 10, 10), # correct state + guessing state
  
  sigma_mu        = as.array(0.4),
  sigma_sigma     = as.array(0.1),
  
  alpha_mu        = c(0.5, 0.25),
  alpha_sigma     = c(0.1, 0.05),
  
  t0_beta         = as.array(5),
  
  init_prob_alpha = rep(5, 2),
  tran_prob_alpha = list(c(8, 2), c(2, 8))
)
#saveRDS(hyperparams, here("saves", "hyperparams.Rds"))

#generated_data <- hmm_later_prior_pred$sample(data = hyperparams, iter_warmup = 0, iter_sampling = n_predictive, fixed_param = TRUE, seed = 2020)
#generated_data$save_object(here("saves", "prior_predictives.Rds"))
generated_data <- readRDS(here("saves", "prior_predictives.Rds"))

# generated_data$summary(variables = c("prop_state"))
# generated_data$summary(variables = c("prop_responses"))
# generated_data$summary(variables = c("mean_rt"))
# 
# mcmc_hist(generated_data$draws(variables = "prop_state"),     facet_args = list(scales = "fixed", dir = "v"))
# mcmc_hist(generated_data$draws(variables = "prop_responses"), facet_args = list(scales = "fixed", dir = "v"))
# mcmc_hist(generated_data$draws(variables = "mean_rt"),        facet_args = list(scales = "fixed", dir = "v"))


# hist(as.vector(generated_data$draws(variables = "rt")[10,1,]), breaks=20)

rt <- generated_data$draws("rt")
responses <- generated_data$draws("responses")
state <- generated_data$draws("state")

# average run length
run_lengths <- sapply(1:dim(state)[1], function(x) mean(rle(as.vector(state[x,,]))$lengths))
summary(run_lengths)
sd(run_lengths)
#

### Parameter recovery ----
parameters <- c("nu_vec[1,1]", "nu_vec[2,1]", "sigma", "alpha", "t0", "init_prob[1]", "tran_prob[1,1]", "tran_prob[2,2]")
generating_parameters <- as_tibble(posterior::as_draws_df(generated_data$draws(parameters)))
hmm_later <- cmdstan_model(stan_file = here("stan", "later", "hmm_later.stan"), include_paths = here())

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

# https://stackoverflow.com/questions/34208564/how-to-hide-or-disable-in-function-printed-message/34208658#34208658
# hide output from Stan
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

### Maximum aposteriori estimation ----
map <- function(stan_data, true_states, max_iter = 50){
  converged <- FALSE
  i <- 1
  while(i <= max_iter && !converged){
    optims <- suppressWarnings(quiet(hmm_later$optimize(data = stan_data, init = initFun, refresh = 0)))
    converged <- length(optims$output_files(include_failed = FALSE)) != 0  
    if(converged){
      # assess label switching
      est_states <- apply(matrix(optims$mle("state_prob"), ncol = 2, byrow = TRUE), 1, which.max)
      p_agree <- mean(est_states == true_states)
      if(p_agree < 0.5) converged <- FALSE
    }
    i <- i + 1
  }
  
  if(converged){
    return(optims$mle(variables = parameters))
  } else{
    nas <- rep(NA, ncol(estimated_parameters))
    names(nas) <- colnames(estimated_parameters)
    return(nas)
  }
}

char2label <- function(x){
  xx <- sprintf("expression(%s)", x)
  eval(parse(text=xx))
}

stan_data <- hyperparams
# stan_data$rt <- as.vector(rt[1,1,])
# stan_data$responses <- as.vector(responses[1,1,])

estimated_parameters <- dplyr::select(generating_parameters, -c(".chain", ".iteration", ".draw"))
parameters_names <- colnames(estimated_parameters)

# fit models (get MAP estimates)
# pb <- dplyr::progress_estimated(n = n_predictive)
# for(i in seq_len(n_predictive)){
#   stan_data$rt <- as.vector(rt[i,1,])
#   stan_data$responses <- as.vector(responses[i,1,])
#   estimated_parameters[i,] <- as.list(map(stan_data, as.vector(state[i,1,])))
#   pb$tick()$print()
# }
# saveRDS(estimated_parameters, file = here("saves", "map_estimates_simulation.Rds"))
estimated_parameters <- readRDS(here("saves", "map_estimates_simulation.Rds"))

parameters_labels <- c("nu[1]^(1)", "nu[1]^(2)", "sigma", "alpha^(1)", "alpha^(2)", "tau", "pi[1]", "rho[11]", "rho[22]")
names(parameters_labels) <- parameters_names

# MAP estimation parameter recovery results:
png(here("figures", "simulations", "parameter_recovery_map.png"), 
    pointsize = 25, width = 750, height = 750)
mean(complete.cases(estimated_parameters))
kk <- sqrt(length(parameters_names))
par(mfrow = c(floor(kk), ceiling(kk)), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0), oma = c(1.5, 1.5, 0.25, 0))
for(p in parameters_names){
  x <- generating_parameters[,p,drop=TRUE]
  y <- estimated_parameters[,p,drop=TRUE]
  lims <- range(c(x, y), na.rm = TRUE)
  plot(x, y, pch = 19, xlab = "", ylab = "", main = "", 
       xlim = lims, ylim = lims, col = adjustcolor("black", alpha = 0.3))
  title(char2label(parameters_labels[p]), line = 1, cex.main = 1.5)
  abline(0, 1, lwd = 3)
  
  corr <- round(cor(x, y, use = "p"), 2)
  text(x = seq(lims[1], lims[2], length.out = 10)[8], 
       y = seq(lims[1], lims[2], length.out = 10)[2], 
       bquote("r" == .(corr))
       )
}
mtext("True", side = 1, outer = TRUE, adj = 0.53, line = -0.5)
mtext("Estimated", side = 2, outer = TRUE, ad = 0.51, line = -0.5)
par(mfrow = c(1, 1))
dev.off()

### Using full Bayes ----
# run run_fit_sbc.sh script

### Simulation based calibration ----
# the distribution of the cumulative probability of true parameter compared to the prior predictive data averaged posterior
# should be uniform

# let's combine all posteriors to get data averaged posterior
saved_fits <- list.files(here("saves", "sbc"))
# we need to order the files
saved_fits <- saved_fits[gsub("sim_", "", saved_fits) %>% gsub(".Rds", "", .) %>% as.integer() %>% order()]

draws <- lapply(seq_along(saved_fits), function(i) {
  file <- here("saves", "sbc", saved_fits[i])
  samples <- readRDS(file)
  rhats <- samples$summary(parameters)$rhat
  labels <- apply(matrix(samples$summary("state_prob")$mean, ncol = 2, byrow = TRUE),1,which.max)
  as_tibble(as_draws_df(samples$draws(parameters))) %>% 
    mutate(.sim = i,
           good_psrf = all(rhats > 0.99) && all(rhats < 1.01),
           good_labels = mean(labels == as.vector(state[i,,])) > 0.5)
})
draws <- bind_rows(draws)

mean(draws$good_labels)
mean(draws$good_psrf)

# for each parameter draw from prior, compute the probability that it is larger than a random draw from a posterior
p_ecdf <- estimated_parameters
for(par in parameters_names){
  for(i in seq_len(n_predictive)){
    p_ecdf[[par]][[i]] <- mean(generating_parameters[[par]][[i]] > subset(draws, .sim == i)[[par]])
  }
}

png(filename = here("figures", "simulations", "sbc_density.png"), pointsize = 25, width = 750, height = 750)
par(mfrow = c(floor(kk), ceiling(kk)), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0), oma = c(0, 1.5, 0.5, 0))
for(par in parameters_names){
  hist(p_ecdf[[par]], breaks = seq(0, 1, 0.1), main = "", xlab = "", ylab = "",
       freq = TRUE, col = "gray")
  title(char2label(parameters_labels[par]), line = 1, cex.main = 1.5)
}
mtext("Frequency", side = 2, outer = TRUE, ad = 0.51, line = -0.5)
par(mfrow = c(1,1))
dev.off()

png(filename = here("figures", "simulations", "sbc_cumulative.png"), pointsize = 25, width = 750, height = 750)
par(mfrow = c(floor(kk), ceiling(kk)), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0), oma = c(0, 1.5, 0.5, 0))
for(par in parameters_names){
  plot(ecdf(p_ecdf[[par]]), main = "", xlab = "", ylab = "")
  title(char2label(parameters_labels[par]), line = 1.2, cex.main = 1.5)
  abline(a = 0, b = 1, col = "red", lty = 2, lwd = 3)
}
mtext("ECDF", side = 2, outer = TRUE, ad = 0.51, line = -0.5)
par(mfrow = c(1,1))
dev.off()

posterior_summaries <- draws %>%
  subset(good_labels) %>%
  pivot_longer(cols = all_of(parameters_names), names_to = "parameter", values_to = "estimated_value") %>%
  group_by(parameter, .sim) %>%
  summarise(mean_post = mean(estimated_value),
            variance_post = var(estimated_value),
            std_dev_post = sd(estimated_value), 
            q10 = quantile(estimated_value, 0.10),
            q25 = quantile(estimated_value, 0.25),
            q75 = quantile(estimated_value, 0.75),
            q90 = quantile(estimated_value, 0.90)) %>%
  ungroup() %>%
  left_join(
    generating_parameters %>% 
      pivot_longer(cols = all_of(parameters_names), names_to = "parameter", values_to = "true_value") %>%
      mutate(.sim = .iteration),
    by = c("parameter", ".sim")
  ) %>%
  select(-".chain", -".iteration", -".draw") %>%
  left_join(
    generating_parameters %>% 
      pivot_longer(cols = all_of(parameters_names), names_to = "parameter", values_to = "true_value") %>%
      group_by(parameter) %>%
      summarise(variance_prior = var(true_value)) %>%
      ungroup(),
    by = "parameter"
  ) %>% 
  mutate(z_score = (mean_post - true_value) / std_dev_post,
         contraction = 1 - variance_post/variance_prior,
         covered50 = true_value > q25 & true_value < q75,
         covered80 = true_value > q10 & true_value < q90)
  
posterior_summaries %>% 
  #subset(z_score < 5 & z_score > -5) %>%
  #subset(contraction > 0 & contraction < 1) %>%
  ggplot(aes(x = contraction, y = z_score)) + 
  geom_point() +
  facet_wrap(.~parameter, scales = "fixed") + 
  theme_bw()

png(here("figures", "simulations", "sensitivity_fixed.png"), pointsize = 25, width = 750, height = 750)
kk <- sqrt(length(parameters_names))
par(mfrow = c(floor(kk), ceiling(kk)), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0), oma = c(1.5, 1.5, 0.25, 0))
xlim <- range(posterior_summaries$contraction)
ylim <- range(posterior_summaries$z_score)
for(p in parameters_names){
  x <- subset(posterior_summaries, subset = parameter == p)[["contraction"]]
  y <- subset(posterior_summaries, subset = parameter == p)[["z_score"]]
  
  plot(x, y, pch = 19, xlab = "", ylab = "", main = "", 
       xlim = xlim, ylim = ylim, col = adjustcolor("black", alpha = 0.3))
  title(char2label(parameters_labels[p]), line = 1, cex.main = 1.5)
  abline(h = 0, v = 0, col = adjustcolor("black", alpha = 0.5), lwd = 2, lty = 2)
}
mtext("Posterior contraction", side = 1, outer = TRUE, adj = 0.53, line = -0.5)
mtext("Posterior z-score", side = 2, outer = TRUE, adj = 0.51, line = -0.5)
par(mfrow = c(1, 1))
dev.off()

# Posterior expectation parameter recovery results:
# quick gg plot
posterior_summaries %>%
  ggplot(aes(x = true_value, y = mean_post)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point() +
  facet_wrap(.~parameter, scales = "free") + 
  theme_bw()

png(here("figures", "simulations", "parameter_recovery_mean.png"), 
    pointsize = 25, width = 750, height = 750)
kk <- sqrt(length(parameters_names))
par(mfrow = c(floor(kk), ceiling(kk)), mar = c(3, 3, 2, 1), mgp = c(2, 1, 0), oma = c(1.5, 1.5, 0.25, 0))
for(p in parameters_names){
  x <- subset(posterior_summaries, subset = parameter == p)[["true_value"]]
  y <- subset(posterior_summaries, subset = parameter == p)[["mean_post"]]
  #browser()
  lims <- range(c(x, y), na.rm = TRUE)
  plot(x, y, pch = 19, xlab = "", ylab = "", main = "", 
       xlim = lims, ylim = lims, col = adjustcolor("black", alpha = 0.3))
  title(char2label(parameters_labels[p]), line = 1, cex.main = 1.5)
  abline(0, 1, lwd = 3)
  
  corr <- round(cor(x, y, use = "p"), 2)
  text(x = seq(lims[1], lims[2], length.out = 10)[8], 
       y = seq(lims[1], lims[2], length.out = 10)[2], 
       bquote("r" == .(corr))
  )
}
mtext("True", side = 1, outer = TRUE, adj = 0.53, line = -0.5)
mtext("Estimated", side = 2, outer = TRUE, adj = 0.51, line = -0.5)
par(mfrow = c(1, 1))
dev.off()


jeffreys <- function(covered, p){
  a <- sum( covered, na.rm = TRUE)
  b <- sum(!covered, na.rm = TRUE)
  
  qbeta(p = p, a + 1/2, b + 1/2)
}
coverage <- posterior_summaries %>% 
  mutate(par_label = parameters_labels[parameter]) %>%
  group_by(par_label) %>%
  summarise(point50 = mean(covered50), lower50 = jeffreys(covered50, 0.025), upper50 = jeffreys(covered50, 0.975),
            point80 = mean(covered80), lower80 = jeffreys(covered80, 0.025), upper80 = jeffreys(covered80, 0.975))
coverage <- coverage[match(parameters_labels, coverage$par_label), ] # reorder rows
# |par_label |         point50| lower50| upper50|         point80| lower80| upper80|
# |:---------|---------------:|-------:|-------:|---------------:|-------:|-------:|
# |nu[1]^(1) |           0.519|   0.487|   0.550|           0.790|   0.764|   0.816|
# |nu[1]^(2) |           0.480|   0.449|   0.512|           0.792|   0.765|   0.817|
# |sigma     |           0.507|   0.475|   0.539|           0.824|   0.799|   0.848|
# |alpha^(1) |           0.487|   0.455|   0.519|           0.783|   0.756|   0.808|
# |alpha^(2) |           0.509|   0.477|   0.541|           0.812|   0.786|   0.836|
# |tau       |           0.501|   0.469|   0.532|           0.814|   0.788|   0.838|
# |pi[1]     |           0.487|   0.455|   0.519|           0.803|   0.777|   0.828|
# |rho[11]   |           0.525|   0.493|   0.557|           0.833|   0.808|   0.856|
# |rho[22]   |           0.507|   0.475|   0.539|           0.799|   0.772|   0.824|

# LaTeX for paper
coverage %>%
  mutate(par_label = sprintf("$\\%s$", par_label)) %>%
  mutate(par_label = gsub("\\[", "_{", par_label)) %>%
  mutate(par_label = gsub("\\]", "}",  par_label)) %>%
  mutate(par_label = gsub("\\(", "{(", par_label)) %>%
  mutate(par_label = gsub("\\)", ")}",  par_label)) %>%
  mutate(cov50 = sprintf("%0.2f [%0.2f, %0.2f]", point50, lower50, upper50),
         cov80 = sprintf("%0.2f [%0.2f, %0.2f]", point80, lower80, upper80)) %>%
  dplyr::select(par_label, cov50, cov80) %>%
  knitr::kable(digits = 2, escape = FALSE, format = "latex", booktabs = TRUE)

#### 
# stan_data <- hyperparams
# stan_data$N_obs <- nrow(tdA)
# stan_data$rt <- tdA$RT/1000
# stan_data$responses <- ifelse(tdA$correct == 1, 1, 2) 
# samples <- hmm_later$sample(data = stan_data, chains = 16, parallel_chains = 4, refresh = 500, init = initFun)
# optims  <- hmm_later$optimize(data = stan_data, init = initFun)
# 
# library(rstan)
# foo <- stan_model(file = here("stan", "later", "hmm_later_prior_pred.stan"))
# dat <- sampling(foo, hyperparams, algorithm = "Fixed_param", chains = 1)
