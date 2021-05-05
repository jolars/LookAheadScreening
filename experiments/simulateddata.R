library(LookAheadScreening)
#library(tibble)
#library(dplyr)
library(tidyr)
library(readr)

printf <- function(...) invisible(cat(sprintf(...)))

g <- expand_grid(
  n = 100,
  p = 20000,
  snr = 3,
  s = 10,
  rho = c(0, 0.4, 0.8),
  screening_type = c("gap_safe", "gap_safe_lookahead"),
  path_length = 100,
  avg_screened = NA,
  time = list(NA),
  screened = list(NA),
  active = list(NA),
  step = list(NA)
)

n_it <- 2

for (i in seq_len(nrow(g))) {
  rho <- g$rho[i]
  screening_type <- g$screening_type[i]
  path_length <- g$path_length[i]

  n <- g$n[i]
  p <- g$p[i]
  snr <- g$snr[i]
  s <- g$s[i]

  avg_screened <- violations <- time <- double(n_it)
  active <- screened <- matrix(NA, nrow = n_it, ncol = path_length)

  printf(
    "%02d/%i n: %4d p: %4d rho: %1.1f %-10s\n",
    i, nrow(g), n, p, rho, screening_type
  )

  for (j in seq_len(n_it)) {
    set.seed(j)

    d <- generateDesign(n, p, family = "gaussian", rho = rho, snr = snr)
    X <- d$X
    y <- d$y

    fit <- lassoPath(
      X,
      y,
      screening_type = screening_type,
      path_length = path_length
    )

    n_lambda <- ncol(fit$beta)

    time[j] <- fit$full_time
    avg_screened[j] <- mean(fit$active / fit$screened)
    screened[j, 1:n_lambda] <- fit$screened
    active[j, 1:n_lambda] <- fit$active

    # # stop if standard error is within +- % of mean
    # if (j > 9) {
    #   time_se <- sd(time[1:j]) / sqrt(j)

    #   if (time_se / mean(time[1:j]) < 0.025) {
    #     break
    #   }
    # }
  }

  time <- time[1:j]
  active <- active[1:j, ]
  screened <- screened[1:j, ]
  avg_screened <- avg_screened[1:j]

  dontuse <- apply(screened, 2, anyNA)

  active <- colMeans(active)
  screened <- colMeans(screened)

  g$n[i] <- n
  g$p[i] <- p
  g$time[i] <- list(time)
  g$avg_screened[i] <- mean(avg_screened)
  g$screened[i] <- list(screened)
  g$active[i] <- list(active)
  g$step[i] <- list(1:path_length)
}

save_rds(g, "results/simulateddata.rds")
