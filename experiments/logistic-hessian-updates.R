library(HessianScreening)

n <- 2000
p <- 10000
rho <- 0.4
dens <- 1

d <-
  generateDesign(n,
                 p,
                 family = "binomial",
                 rho = rho,
                 density = dens)

X <- d$X
y <- d$y

fit_auto <-
  lassoPath(
    X,
    y,
    family = "binomial",
    screening_type = "hessian",
    log_hessian_update_type = "auto",
    log_hessian_auto_threshold = ceiling(dens*min(n, p)),
    verbosity = 0
  )
fit_full <-
  lassoPath(
    X,
    y,
    family = "binomial",
    screening_type = "hessian",
    log_hessian_update_type = "full"
  )
fit_approx <-
  lassoPath(
    X,
    y,
    family = "binomial",
    screening_type = "hessian",
    log_hessian_update_type = "approx"
  )
fit_working <-
  lassoPath(X,
            y,
            family = "binomial",
            screening_type = "working")

fit_auto$full_time
fit_full$full_time
fit_approx$full_time
fit_working$full_time

library(tikzDevice)
tikz("log-hessian.pdf")
par(mfrow = c(2, 2))

plot(
  fit_auto$passes,
  type = "l",
  ylab = "passes",
  ylim = c(
    0,
    max(
      fit_full$passes,
      fit_approx$passes,
      fit_auto$passes,
      fit_working$passes
    )
  )
)
lines(fit_full$passes, col = "dark orange")
lines(fit_approx$passes, col = "aquamarine3")
lines(fit_working$passes, col = "magenta")
legend(
  "topleft",
  legend = c("auto", "full", "approx", "working"),
  lwd = c(1, 1, 1, 1),
  col = c("black", "dark orange", "aquamarine3", "magenta")
)

plot(
  fit_auto$it_time,
  type = "l",
  ylab = "iteration time",
  ylim = c(0, max(
    c(fit_auto$it_time, fit_full$it_time, fit_approx$it_time)
  ))
)
lines(fit_full$it_time, col = "dark orange")
lines(fit_approx$it_time, col = "aquamarine3")
lines(fit_working$it_time, col = "magenta")
legend(
  "topleft",
  legend = c("auto", "full", "approx", "working"),
  lwd = c(1, 1, 1, 1),
  col = c("black", "dark orange", "aquamarine3", "magenta")
)

plot(
  fit_auto$active,
  type = "l",
  ylab = "screened",
  lty = 2,
  ylim = c(0, max(
    c(
      fit_auto$screened,
      fit_full$screened,
      fit_approx$screened,
      fit_working$screened
    )
  ))
)
lines(fit_auto$screened)
lines(fit_full$screened, col = "dark orange")
lines(fit_approx$screened, col = "aquamarine3")

legend(
  "topleft",
  legend = c("auto", "full", "approx"),
  lwd = c(1, 1, 1),
  col = c("black", "dark orange", "aquamarine3")
)

plot(
  cumsum(fit_auto$it_time),
  type = "l",
  ylab = "cumulative time",
  ylim = c(0, max(c(
    cumsum(fit_auto$it_time),
    cumsum(fit_full$it_time),
    cumsum(fit_approx$it_time),
    cumsum(fit_working$it_time)
  )))
)
lines(cumsum(fit_full$it_time), col = "dark orange")
lines(cumsum(fit_approx$it_time), col = "aquamarine3")
lines(cumsum(fit_working$it_time), col = "magenta")
legend(
  "topleft",
  legend = c("auto", "full", "approx", "working"),
  lwd = c(1, 1, 1, 1),
  col = c("black", "dark orange", "aquamarine3", "magenta")
)
mtext(
  paste0("n = ", n, ", p = ", p, ", rho = ", rho, ", density = ", dens),
  line = -2,
  outer = TRUE
)
par(mfrow = c(1, 1))
dev.off()
