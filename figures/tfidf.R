
library(ggplot2)
library(dplyr)
library(tidyr)
library(tikzDevice)

source("R/utils.R")

d <- readRDS("results/simulateddata.rds")

theme_set(theme_minimal(base_size = 9))

cols <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7"
)

options(
  tikzDocumentDeclaration =
    "\\documentclass[10pt]{article}\n\\usepackage{newtxtext,newtxmath}\n"
)

d2 <-
  d %>%
  mutate(
    np = paste0("$n=", n, "$, $p=", p, "$"),
    screening_type = recode(
      screening_type,
      "hessian" = "Hessian",
      "working" = "Working",
      "gap_safe" = "Gap Safe",
      "edpp" = "EDPP",
      "strong" = "Strong"
    )
  ) %>%
  filter(screening_type != "Working") %>%
  drop_na() %>%
  unnest(c(screened, active, step)) %>%
  group_by(family, rho, scenario, screening_type) %>%
  mutate(screened = screened / p, active = active / p)

rho_labeller <- function(labels, multi_line = TRUE, sep = ":", ...) {
  value <- label_value(labels, multi_line = multi_line)

  out <- paste0("$\\rho = ", value, "$")
  #  list(unname(unlist(out)))
  out
}

d2_gaussian <- filter(d2, family == "gaussian")
d2_gaussian_active <-
  d2_gaussian %>%
  ungroup() %>%
  group_by(family, rho, np, step) %>%
  summarize(avg_active = mean(active, na.rm = TRUE))

f1 <- "figures/simulateddata-efficiency-gaussian.tex"
tikz(f1, width = 5.6, height = 2.5, standAlone = TRUE)
ggplot(d2_gaussian, aes(step)) +
  facet_grid(np ~ rho,
    labeller = labeller(
      np = label_value,
      rho = rho_labeller
    )
  ) +
  geom_line(aes(x = step, y = avg_active),
    linetype = 2,
    data = d2_gaussian_active
  ) +
  geom_line(aes(y = screened, col = screening_type, group = screening_type)) +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = cols) +
  labs(y = "Fraction of Predictors", x = "Step")
dev.off()
renderPdf(f1)

d2_binomial <- filter(d2, family == "binomial")
d2_binomial_active <-
  d2_binomial %>%
  ungroup() %>%
  group_by(family, rho, np, step) %>%
  summarize(avg_active = mean(active, na.rm = TRUE))

f2 <- "figures/simulateddata-efficiency-binomial.tex"
tikz(f2, width = 5.6, height = 2.5, standAlone = TRUE)
ggplot(d2_binomial, aes(step)) +
  facet_grid(np ~ rho,
    labeller = labeller(
      np = label_value,
      rho = rho_labeller
    )
  ) +
  geom_line(aes(x = step, y = avg_active),
    linetype = 2,
    data = d2_binomial_active
  ) +
  geom_line(aes(y = screened, col = screening_type, group = screening_type)) +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = cols[-1]) +
  labs(y = "Fraction of Predictors", x = "Step")
dev.off()
renderPdf(f2)
