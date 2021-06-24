library(tibble)
library(tidyr)
library(dplyr)
library(forcats)
library(tikzDevice)
library(ggplot2)
library(readr)
library(here)

source(here("R", "utils.R"))

theme_set(theme_minimal(base_size = 9))

d1 <- read_rds(here("results", "simulateddata.rds")) %>%
  select(-(screened:step)) %>%
  unnest(time) %>%
  mutate(
    screening_type = recode(
      screening_type,
      "gap_safe" = "Standard",
      "gap_safe_lookahead" = "Look-Ahead"
    ),
    snr = as.factor(snr)
  )

cols <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7"
)

file <- here("paper", "figures", "simulated-data-timings.tex")
tikz(file, width = 4.5, height = 1.8, standAlone = TRUE)
ggplot(d1, aes(
  snr,
  time,
  fill = screening_type
)) +
  geom_boxplot(outlier.alpha = 0.7) +
  scale_fill_manual(values = cols) +
  labs(fill = NULL, x = "Signal-to-Noise Ratio", y = "Time (s)")
dev.off()

renderPdf(file)
