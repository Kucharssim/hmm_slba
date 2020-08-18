library(here)
source(here("R", "expose_stan_functions.R"))

png(here("figures", "dlater_single.png"), width = 800, height = 500, pointsize = 20)
curve(dlater_single(rt=x, 1, 0.5, 0.5, 0.5), from = 0, to = 2, lwd = 3, n = 1000,
      bty = "l", xlab = "Time (sec)", ylab = "Density")
dev.off()


curve(dlater(rt=x, 1, c(0.6, 0.4), c(0.2, 0.2), c(0.5, 0.5), 0.1), from = 0, to = 3)
