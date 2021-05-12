
library(LookAheadScreening)

# d <- readRDS(file.path("data", paste0("leukemia-train", ".rds")))
d <- generateDesign(100, 10, s = 1, density = 1)
X <- d$X
y <- d$y
n <- nrow(X)
p <- ncol(X)
tol_gap <- 1e-5
tol_infeas <- 1e-4
verbosity <- 0

fit <- lassoPath(
    X,
    y,
    screening_type = "gap_safe",
    standardize = TRUE,
    verbosity = 1,
    tol_gap = tol_gap,
    tol_infeas = tol_infeas,
    path_length = 100
)

fit_lookahead <- lassoPath(
    X,
    y,
    screening_type = "gap_safe_lookahead",
    standardize = TRUE,
    verbosity = 1,
    tol_gap = tol_gap,
    tol_infeas = tol_infeas,
    path_length = 100
)

library(tibble)
library(ggplot2)

d <- tibble(
    time = c(fit_lookahead$cd_time, fit$cd_time),
    passes = c(fit_lookahead$passes, fit$passes),
    screened = c(fit_lookahead$screened, fit$screened),
    step = c(seq_along(fit_lookahead$cd_time), seq_along(fit$cd_time)),
    type = rep(c("lookahead", "standard"),
        times = c(length(fit_lookahead$cd_time), length(fit$cd_time))
    )
)

ggplot(d, aes(step, screened, col = type)) +
    geom_line()

fit$full_time
fit_lookahead$full_time

fit_lookahead$lookahead
fit_lookahead$violations
fit$violations

any(fit_lookahead$violations > 0)

fit_lookahead$lookahead[1:30, 1:39]
