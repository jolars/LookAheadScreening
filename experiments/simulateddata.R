
library(HessianScreening)
library(RColorBrewer)
library(tibble)
library(dplyr)
library(tidyr)
library(progress)

g <- expand_grid(
  np = list(c(5e4, 1e2), c(1e2, 5e4)),
  n = NA,
  p = NA,
  rho = c(0, 0.3, 0.6),
  family = c("binomial", "gaussian"),
  screening_type = c("hessian", "working", "gap_safe", "edpp"),
  path_length = 100,
  avg_screened = NA,
  avg_violations = NA,
  time = list(NA),
  screened = list(NA),
  active = list(NA)
)

n_it <- 20

pb <- progress_bar$new(
  format = "  simulating [:bar] :percent eta: :eta",
  total = nrow(g) * n_it,
  clear = FALSE,
  width = 80
)

for (i in seq_len(nrow(g))) {
  np <- g$np[i]
  n <- np[[1]][1]
  p <- np[[1]][2]
  rho <- g$rho[i]
  family <- g$family[i]
  screening_type <- g$screening_type[i]
  path_length <- g$path_length[i]

  if (family == "binomial" && screening_type == "edpp") {
    next
  }

  avg_screened <- violations <- time <- double(n_it)
  active <- screened <- matrix(NA, nrow = n_it, ncol = path_length)

  for (j in seq_len(n_it)) {
    set.seed(j)
    pb$tick()

    d <- generateDesign(n, p, family = family, rho = rho)
    X <- d$X
    y <- d$y

    fit <- lassoPath(
      X,
      y,
      family = family,
      screening_type = screening_type,
      path_length = path_length
    )

    n_lambda <- length(fit$lambda)

    time[j] <- fit$full_time
    avg_screened[j] <- mean(fit$active / fit$screened)
    violations[j] <- sum(fit$violations)
    screened[j, 1:n_lambda] <- fit$screened
    active[j, 1:n_lambda] <- fit$active
  }

  keep <- !apply(screened, 2, anyNA)

  active <- active[, keep]
  screened <- screened[, keep]

  g$n[i] <- n
  g$p[i] <- p
  g$family[i] <- family
  g$time[i] <- list(time)
  g$avg_screened[i] <- mean(avg_screened)
  g$avg_violations[i] <- mean(violations)
  g$screened[i] <- list(screened)
  g$active[i] <- list(active)
}

saveRDS(g, "results/simulateddata.rds")
