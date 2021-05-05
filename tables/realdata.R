library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(forcats)
library(stringr)

d_raw <- readRDS("results/realdata.rds")

d <-
  d_raw %>%
  drop_na(time) %>%
  mutate(
    screening_type = recode(
      screening_type,
      "working" = "Working",
      "hessian" = "Hessian",
      "gap_safe" = "GapSafe",
      "edpp" = "EDPP"
    ),
    family = recode(
      family,
      "gaussian" = "Gaussian",
      "binomial" = "Binomial"
    ),
    dataset = str_remove(dataset, "(-train|-test)"),
    dataset = recode(
      dataset,
      "YearPredictionMSD" = "YearPredMSD"
    )
  ) %>%
  select(dataset, family, n, p, density, screening = screening_type, time) %>%
  arrange(family, dataset, screening) %>%
  pivot_wider(names_from = "screening", values_from = "time")

filter(d, family == "Gaussian") %>%
  select(-family) %>%
  write_csv("tables/realdata-timings-gaussian.csv")

filter(d, family == "Binomial") %>%
  select(-family) %>%
  write_csv("tables/realdata-timings-binomial.csv")
