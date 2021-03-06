library(LookAheadScreening)
library(tidyr)
library(readr)

printf <- function(...) invisible(cat(sprintf(...)))

g <- expand_grid(
  n = 100,
  p = 50000,
  snr = c(0.1, 1, 6),
  s = c(5),
  screening_type = c("gap_safe", "gap_safe_lookahead"),
  path_length = c(100),
  avg_screened = NA,
  time = list(NA),
  screened = list(NA),
  active = list(NA),
  step = list(NA)
)

n_it <- 50

for (i in seq_len(nrow(g))) {
  screening_type <- g$screening_type[i]
  path_length <- g$path_length[i]

  n <- g$n[i]
  p <- g$p[i]
  snr <- g$snr[i]
  s <- g$s[i]

  avg_screened <- violations <- time <- double(n_it)
  active <- screened <- matrix(NA, nrow = n_it, ncol = path_length)

  printf(
    "%02d/%i n: %4d p: %4d snr: %.1f %-10s\n",
    i, nrow(g), n, p, snr, screening_type
  )

  for (j in seq_len(n_it)) {
    set.seed(j)

    d <- generateDesign(n, p, snr = snr)
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
  }

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

# library(dplyr)

# unnest(g, time) %>%
#   group_by(snr, path_length, screening_type) %>%
#   summarize(mean_time = mean(time))

write_rds(g, "results/simulateddata.rds")
