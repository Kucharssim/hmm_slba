library(here)
# plot LBA ----

meandrift1 <- 0.6
meandrift2 <- 0.4
sddrift <- 0.16
startpoint <- 0.15
boundary <- 0.37


png(here("figures", "lba.png"), width = 600, height = 500)
par(oma = c(0,2, 0, 0), cex = 1.5)
plot(0, type = "n", bty = "n", xlim = c(-0.02, 0.8), ylim = c(0, 0.5), 
     axes = FALSE, xlab = "Decision Time", ylab = "", yaxs = "i", xaxs = "i")
rect(-0.02, 0, 0, startpoint, col = "gray")
axis(1, lwd = 2)
axis(2, las = 1, lwd = 2,
     at = c(0, startpoint, boundary, 0.5), 
     labels = c(0, "Start point \n upper bound", "Decision \n boundary", ""))
abline(h = boundary, lty = 2, lwd = 1.5)
start <- 0.3*startpoint
arrows(x0 = 0, y0 = start, x1 = (boundary-start)/meandrift1, y = boundary, lwd = 3, length = 0.15)
text(0.035, start, "a")
text(0.2, 0.13, "b")
text((boundary-start)/meandrift1, boundary+0.025, "c")
dev.off()