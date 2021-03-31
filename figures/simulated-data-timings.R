library(tibble)
library(tidyr)
library(dplyr)
library(forcats)
library(lattice)
library(latticeExtra)
library(tactile)
library(tikzDevice)
library(ggplot2)

thm <- tactile.theme(c(12, 8),
  plot.symbol = list(col = 1),
  superpose.line = list(lty = 1:3)
)

d_raw <- readRDS("results/simulateddata.R") %>%
  mutate(
    screening_type = recode(
      screening_type,
      "hessian" = "Hessian",
      "working" = "Working",
      "edpp" = "EDPP",
      "gap_safe" = "Gap-Safe"
    ),
    rho = as.factor(rho),
    np = paste(n, p, sep = " x ")
  ) %>%
  select(np, n, p, rho, family, screening_type, time) %>%
  unnest(time)

d1 <-
  d_raw %>%
  mutate(
    screening_type = as.factor(screening_type),
    family = as.factor(family)
  ) %>%
  group_by(np, rho, family, screening_type) %>%
  summarize(
    meantime = mean(time),
    se = sd(time / n())
  ) %>%
  mutate(
    hi = meantime + se,
    lo = meantime - se
  )

d2_gaussian <-
  filter(d1, family == "gaussian")

library(qualpalr)

cols <- qualpal(4)

cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tikz("figures/simulateddata-gaussian-timings.tex", width = 5, height = 3)
xyplot(
  meantime ~ rho | np,
  groups = screening_type,
  data = d2_gaussian,
  type = "b",
  auto.key = list(lines = TRUE, points = FALSE),
  prepanel = prepanel.ci,
  lower = d2_gaussian$lo,
  upper = d2_gaussian$hi,
  panel = function(...) {
    panel.ci(..., alpha = 0.2, grid = TRUE)
    panel.xyplot(...)
  }
)
dev.off()

library(ggthemes)

tikz("figures/simulateddata-gaussian-timings.tex", width = 5, height = 3)
pl <- ggplot(d2_gaussian, aes(rho,
  meantime,
  ymin = lo,
  ymax = hi,
  group = screening_type,
  color = screening_type,
  fill = screening_type
)) +
  geom_col(position = "dodge", col = 1) +
  facet_wrap("np") +
  theme_minimal() +
  # scale_fill_manual(values = cols) +
  # scale_fill_colorblind() +
  scale_fill_tableau() +
  labs(color = "Screening", fill = "Screening")
dev.off()

d_gaussian <- filter(d_raw, family == "gaussian", p == 10000)

# tikz("figures/simulateddata-gaussian-timings.tex", width = 5, height = 3)
bwplot2(
  rho ~ time,
  col = 1,
  grid = TRUE,
  groups = screening_type,
  data = d_gaussian,
  auto.key = list(space = "right"),
  ylab = expression(rho),
  xlab = "Time",
  par.settings = thm
)
# dev.off()