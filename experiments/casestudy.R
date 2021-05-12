
library(LookAheadScreening)
library(tidyr)
library(tibble)
library(dplyr)
library(readr)

d <- read_rds("data/leukemia-train.rds")

X <- d$X
y <- d$y

n <- nrow(X)
p <- ncol(X)

fit <- lassoPath(
  X,
  y,
  screening_type = "gap_safe_lookahead",
  verbosity = 0
)

lookahead <- fit$lookahead

d <- matrix(0, p, ncol(fit$lookahead))

for (j in 1:p) {
  d[j, 1:(lookahead[j, 1] + 1)] <- TRUE
}

colnames(d) <- 1:ncol(d)

d1 <- tibble(
  step = seq_len(ncol(lookahead)),
  screened = colSums(d) / p
)

write_rds(d1, "results/casestudy-summary.rds")

q <- 20

set.seed(852)

d2 <-
  as.data.frame(d) %>%
  as_tibble() %>%
  sample_n(q) %>%
  rowid_to_column() %>%
  pivot_longer(-rowid) %>%
  mutate(
    name = as.integer(name),
    value = factor(value, labels = c("Included", "Discarded"))
  )

write_rds(d2, "results/casestudy-full.rds")
