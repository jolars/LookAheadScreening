library(LookAheadScreening)
library(tidyr)
library(tibble)
library(Matrix)

printf <- function(...) invisible(cat(sprintf(...)))

datasets <- c(
  "arcene",
  "colon-cancer",
  "e2006-tfidf-train",
  # "e2006-log1p-train",
  # "rcv1-train",
  # "news20"
)

g <- expand_grid(
  dataset = datasets,
  screening_type = c("gap_safe", "gap_safe_lookahead"),
  n = NA,
  p = NA,
  density = NA,
  time = NA,
  avg_screened = NA,
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

  printf("%02d/%i %-10.10s %s\n", i, nrow(g), g$dataset[i], screening_type)

  n_it <- 5

  time <- double(n_it)

  for (k in seq_len(n_it)) {
    fit <- lassoPath(
      X,
      y,
      screening_type = screening_type,
      verbosity = 0
    )

    time[k] <- fit$full_time
  }

  g$n[i] <- n
  g$p[i] <- p
  g$time[i] <- mean(time)
  g$density[i] <- dens
  g$screened[i] <- list(fit$screened)
  g$active[i] <- list(fit$active)
}

cat("DONE!\n")

saveRDS(g, "results/realdata.rds")
