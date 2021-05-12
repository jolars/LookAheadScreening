
library(readr)
library(ggplot2)
library(tikzDevice)

source("R/utils.R")

theme_set(theme_minimal(base_size = 9))

d <- read_rds("results/casestudy-full.rds")

pl <-
  ggplot(d2, aes(name, rowid, fill = value)) +
  geom_tile(show.legend = FALSE, col = "white") +
  scale_fill_manual(values = c("grey", "steelblue4")) +
  coord_fixed() +
  scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 25)) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90), expand = expansion(add = 1)) +
  labs(x = "Step", y = "Predictor")

file <- "paper/figures/casestudy.tex"

tikz(file, width = 4.7, height = 4, standAlone = TRUE)
pl
dev.off()
renderPdf(file)

d2 <- read_rds("results/casestudy-summary.rds")
