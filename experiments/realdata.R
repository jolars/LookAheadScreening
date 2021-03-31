library(HessianScreening)
library(tidyr)
library(tibble)
library(Matrix)

printf <- function(...) invisible(cat(sprintf(...)))

datasets <- c(
  "arcene",
  "cadata",
  "dorothea",
  "gisette-train",
  "colon-cancer",
  "leukemia-train",
  "ijcnn1-train",
  "YearPredictionMSD-train",
  "madelon-train",
  "e2006-tfidf-train",
  "e2006-log1p-train",
  "rcv1-train",
  "news20"
)

g <- expand_grid(
  dataset = datasets,
  screening_type = c("working", "hessian", "gap_safe", "edpp"),
  family = NA,
  n = NA,
  p = NA,
  density = NA,
  time = NA,
  total_violations = NA,
  avg_screened = NA,
  violations = list(NA),
  screened = list(NA),
  active = list(NA)
)

for (i in seq_len(nrow(g))) {
  d <- readRDS(file.path("data", paste0(g$dataset[i], ".rds")))
  screening_type <- g$screening_type[i]

  X <- d$X
  y <- d$y

  n <- nrow(X)
  p <- ncol(X)

  dens <- ifelse(inherits(X, "sparseMatrix"), Matrix::nnzero(X) / length(X), 1)

  family <- if (length(unique(d$y)) == 2) "binomial" else "gaussian"

  if (family == "binomial" && screening_type == "edpp") {
    next
  }

  log_hessian_update_type <-
    ifelse((1 - dens) * min(n, p) / max(n, p) * 10 < 0.1, "full", "approx")

  printf("%02d/%i %-10.10s %s\n", i, nrow(g), g$dataset[i], screening_type)

  n_it <- 1

  time <- double(n_it)

  for (k in seq_len(n_it)) {
    fit <- lassoPath(
      X,
      y,
      family = family,
      screening_type = screening_type,
      verbosity = 0,
      log_hessian_update_type = log_hessian_update_type
    )

    time[k] <- fit$full_time
  }

  g$n[i] <- n
  g$p[i] <- p
  g$family[i] <- family
  g$time[i] <- mean(time)
  g$density[i] <- dens
  g$total_violations[i] <- sum(fit$violations)
  g$avg_screened[i] <- mean(fit$active / fit$screened)
  g$violations[i] <- list(fit$violations)
  g$screened[i] <- list(fit$screened)
  g$active[i] <- list(fit$active)
}

cat("DONE!\n")

saveRDS(g, "results/realdata.rds")
