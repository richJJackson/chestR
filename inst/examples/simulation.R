# Self-contained simulation example for chestR
# Run with: source(system.file("examples", "simulation.R", package = "chestR"))

library(survival)

if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("Install mvtnorm to run this example: install.packages('mvtnorm')")
}

set.seed(181125)
N <- 400

a1 <- rbinom(N, 1, 0.3)
a2 <- rbinom(N, 1, 0.6)
trt <- rbinom(N, 1, 0.5)

sig <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
b <- mvtnorm::rmvnorm(N, c(0, 0), sig)

centre <- c(0.5, 0.5)
centre.sd <- matrix(c(0.2, 0, 0, 0.2), 2, 2)
biom.lp <- mvtnorm::dmvnorm(b, centre, centre.sd)
trt.biom <- -biom.lp * trt

lp <- 0.05 + log(0.5) * a1 + log(0.7) * a2 + trt.biom
st <- rexp(N, exp(lp))
cen <- as.numeric(st < 4)
st[cen == 0] <- 4

base <- coxph(Surv(st, cen) ~ a1 + a2 + trt)
cr <- chestr(base, b, method = "legacy", kern.adj = 2, grid.size = 50)

nev <- sum(cen)
cr$propInf <- cr$eff.event / nev
cex.w <- 4 * cr$propInf / max(cr$propInf)

message("Example complete: ", nrow(cr), " grid points evaluated.")

if (interactive()) {
  plot(b[, 1], b[, 2], xlab = "biomarker 1", ylab = "biomarker 2",
       main = "Local treatment effect (simulation)")
  points(cr$biom1, cr$biom2, cex = cex.w, pch = 20, col = grDevices::heat.colors(10)[cut(cr$trt, 10)])
  points(b[, 1], b[, 2], pch = ".", cex = 3)
}
