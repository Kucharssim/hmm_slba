library(cmdstanr)
set_cmdstan_path(readRDS("path_to_cmdstan.Rds"))
cmdstan_version()
library(here)

### Parameter recovery ----
# we are interested in the following 9 parameters:
parameters <- c("nu_vec[1,1]", "nu_vec[2,1]", "sigma", "alpha", "t0", "init_prob[1]", "tran_prob[1,1]", "tran_prob[2,2]")
generating_parameters <- as_tibble(posterior::as_draws_df(generated_data$draws(parameters)))
# load the stan model to fit the data
hmm_later <- cmdstan_model(stan_file = here("stan", "later", "hmm_later.stan"), include_paths = here())
# load the hyperparameters (for prior specs)
hyperparams <- readRDS(here("saves", "hyperparams.Rds"))

# function that generates random starting values using the priors -----
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

### Maximum aposteriori estimation function ----
# 1) fit model
# 2) check whether states were switched, if yes, refit the model
# 3) if model does not converge after 50 repetitions, stop fitting and return NA for each parameter
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

# prepare a list of hyperparameters + data to be passed to the stan model
stan_data <- hyperparams
stan_data$rt <- NA
stan_data$responses <- NA

# prepare a data frame to store parameter estimates
estimated_parameters <- dplyr::select(generating_parameters, -c(".chain", ".iteration", ".draw"))
parameters_names <- colnames(estimated_parameters)

# fit models (get MAP estimates): uncomment to redo the simulation
# pb <- dplyr::progress_estimated(n = n_predictive)
# for(i in seq_len(n_predictive)){
#   stan_data$rt <- as.vector(rt[i,1,])
#   stan_data$responses <- as.vector(responses[i,1,])
#   estimated_parameters[i,] <- as.list(map(stan_data, as.vector(state[i,1,])))
#   pb$tick()$print()
# }
# saveRDS(estimated_parameters, file = here("saves", "map_estimates_simulation.Rds"))
estimated_parameters <- readRDS(here("saves", "map_estimates_simulation.Rds"))


### plot results ----
char2label <- function(x){
  xx <- sprintf("expression(%s)", x)
  eval(parse(text=xx))
}
parameters_labels <- c("nu[1]^(1)", "nu[1]^(2)", "sigma", "alpha^(1)", "alpha^(2)", "tau", "pi[1]", "rho[11]", "rho[22]")
names(parameters_labels) <- parameters_names

# MAP estimation parameter recovery results ----
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