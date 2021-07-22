library(here)
source(here("R", "expose_stan_functions.R"))

nu_1 <- c(0.6, 0.4)
nu_2 <- c(0.51, 0.49)
alpha_1 <- c(0.35, 0.35)
alpha_2 <- c(0.18, 0.18)
sigma <- c(0.2,0.2)
tau <- 0.2

pi_1 <- 0.46
rho_11 <- 0.922
rho_22 <- 0.89

png(here("figures", "slba-hmm.png"), width = 5, height = 8, res = 300, units = "in")

layout.matrix <- matrix(c(1, 3, 2, 4, 5, 5, 6, 6), nrow = 4, byrow = TRUE)
layout(layout.matrix, heights = c(1, 3, 0.5, 1), widths = c(1, 1))
par(oma = c(2, 2, 1, 2), mar = c(0, 1, 0, 1), cex = 1.2)

# controlled state
df <- data.frame(correct = rlater_single(1000, nu_1[1], sigma[1], alpha_1[1], tau), 
                 error   = rlater_single(1000, nu_1[2], sigma[2], alpha_1[2], tau))
df$rt       <- apply(df[,1:2], 1, min)
df$response <- apply(df[,1:2], 1, which.min)
plot(0, type = "n", xlim = c(0, 1.5), ylim = c(0, 3), axes = FALSE, xlab = "", ylab = "")
axis(2, at = seq(0, 3, by = 1))
mtext("Density", side = 2, cex = 1.5, line = 2)
curve(dlater(rt=x, response = 1, nu = nu_1, sigma = sigma, alpha = alpha_1, t0 = tau), 
      from = 0, to = 1.5, col = "steelblue", lwd = 3, add = TRUE)
curve(dlater(rt=x, response = 2, nu = nu_1, sigma = sigma, alpha = alpha_1, t0 = tau), 
      from = 0, to = 1.5, col = "firebrick", lwd = 3, add = TRUE)
rug(df$rt[df$response==1], col = "steelblue", ticksize = 0.01)
rug(df$rt[df$response==2], col = "firebrick", ticksize = 0.01)



plot(0, type = "n", bty = "n",
     xlim = c(0, 1.4), ylim = c(0, 0.45), axes = FALSE, xlab = "", ylab = "")
axis(1, at = seq(0, 1.5, by = 0.5))
axis(2, at = seq(0, 0.4, by = 0.1))
mtext("Evidence", side = 2, line = 2, cex = 1.5)

# alpha
abline(h = alpha_1, lwd = 2.5, lty = 2)
text(1.3, alpha_1+0.02, expression(alpha^(1)))

# drifts
for(i in 1:nrow(df)) {
  segments(tau, 0, df$correct[i], alpha_1[1], 
           col = adjustcolor("steelblue", 0.02))
  segments(tau, 0, df$error[i], alpha_1[2],
           col = adjustcolor("firebrick", 0.02))
}

arrows(tau, 0, alpha_1[1]/nu_1[1]+tau, alpha_1[1], lwd = 2.5, col = "steelblue", length = 0.1)
arrows(tau, 0, alpha_1[2]/nu_1[2]+tau, alpha_1[2], lwd = 2.5, col = "firebrick", length = 0.1)
text(0.4, 0.25, expression(nu[1]^(1)))
text(0.6, 0.08,  expression(nu[2]^(1)))

# tau
arrows(0, 0, tau, 0, code = 3, angle = 90, length = 0.05, lwd = 1.5)
text(tau/2, 0.02, expression(tau))

# guessing state
df <- data.frame(correct = rlater_single(1000, nu_2[1], sigma[1], alpha_2[1], tau), 
                 error   = rlater_single(1000, nu_2[2], sigma[2], alpha_2[2], tau))
df$rt       <- apply(df[,1:2], 1, min)
df$response <- apply(df[,1:2], 1, which.min)

plot(0, type = "n", xlim = c(0, 1.5), ylim = c(0, 3), axes = FALSE, xlab = "", ylab = "")
curve(dlater(rt=x, response = 1, nu = nu_2, sigma = sigma, alpha = alpha_2, t0 = tau), 
      from = 0, to = 1.5, col = "steelblue", lwd = 3, add = TRUE)
curve(dlater(rt=x, response = 2, nu = nu_2, sigma = sigma, alpha = alpha_2, t0 = tau), 
      from = 0, to = 1.5, col = "firebrick", lwd = 3, add = TRUE)
rug(df$rt[df$response==1], col = "steelblue", ticksize = 0.01)
rug(df$rt[df$response==2], col = "firebrick", ticksize = 0.01)



plot(0, type = "n", bty = "n",
     xlim = c(0, 1.4), ylim = c(0, 0.5), axes = FALSE, xlab = "", ylab = "")
axis(1, at = seq(0, 1.5, by = 0.5))

# alpha
abline(h = alpha_2, lwd = 2.5, lty = 2)
text(1.3, alpha_2+0.02, expression(alpha^(2)))

# drifts
for(i in 1:nrow(df)) {
  segments(tau, 0, df$correct[i], alpha_2[1], 
           col = adjustcolor("steelblue", 0.02))
  segments(tau, 0, df$error[i], alpha_2[2],
           col = adjustcolor("firebrick", 0.02))
}

arrows(tau, 0, alpha_2[1]/nu_2[1]+tau, alpha_2[1], lwd = 2.5, col = "steelblue", length = 0.1)
arrows(tau, 0, alpha_2[2]/nu_2[2]+tau, alpha_2[2], lwd = 2.5, col = "firebrick", length = 0.1)
text(0.2, 0.1,  expression(nu[1]^(2)))
text(0.5, 0.05, expression(nu[2]^(2)))

# tau
arrows(0, 0, tau, 0, code = 3, angle = 90, length = 0.05, lwd = 1.5)
text(tau/2, 0.02, expression(tau))

legend(x = 0, y = 0.4, 
       legend = c("Correct", "Incorrect"), 
       col = c("steelblue", "firebrick"), 
       lty = c(1, 1), lwd = c(5, 5), cex = 0.7, box.lty = 0)
# x-axis plot
plot(0, type = "n", axes = FALSE, xlim = c(-0.25, 0.25), ylim = c(-0.05, 0.4))
text(0, 0, "Response time")

# state plot
plot(0, type = "n", axes = FALSE, xlim = c(-0.25, 1.25), ylim = c(-1, 0.5))
text(0, 0, "Controlled \n state")
arrows(0, -0.6, 0, -0.3, length = 0.1, lwd = 2.5)
text(0, -0.8, expression(pi[1]))
text(1, 0, "Guessing \n state")
arrows(1, -0.6, 1, -0.3, length = 0.1, lwd = 2.5)
text(1, -0.8, expression(pi[2]))

arrows(0.25, 0.1, 0.75, 0.1, length = 0.1, lwd = 2.5)
text(0.5, 0.25, expression(rho[12]))
arrows(0.75, -0.1, 0.25, -0.1, length = 0.1, lwd = 2.5)
text(0.5, -0.25, expression(rho[21]))

dev.off()