library(ggplot2)
library(dplyr)
library(tidyr)
library(lemon)
library(tikzDevice)
library(knitr)

source("R/utils.R")

theme_set(theme_minimal(base_size = 9))

d <- readRDS("results/efficiency-realdata.rds")

# https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend2 <- function(p) {
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(
    facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]])
  )
  empty.facet.panels <- facet.panels[empty.facet.panels]
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  reposition_legend(p, "center", panel = names)
}

d2 <-
  d %>%
  mutate(screening_type = recode(screening_type,
    "edpp" = "EDPP",
    "gap_safe" = "Gap Safe",
    "hessian" = "Hessian",
    "strong" = "Strong"
  )) %>%
  drop_na() %>%
  unnest(c(screened, active)) %>%
  group_by(dataset, screening_type) %>%
  mutate(screened = screened / p, active = active / p) %>%
  mutate(step = seq_along(screened))

cols <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7"
)

pl <- ggplot(d2, aes(step)) +
  facet_wrap(~dataset, nrow = 2) +
  geom_line(aes(y = active), linetype = 2) +
  geom_line(aes(y = screened, col = screening_type)) +
  scale_y_continuous(
    oob = function(x, ...) {
      x
    },
    limits = c(0, 0.06)
  ) +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = cols) +
  labs(x = "Step", y = "Active Predictors")

file <- "figures/realdata-efficiency.tex"

tikz(file, width = 5.6, height = 2.5, standAlone = TRUE)
shift_legend2(pl)
dev.off()

renderPdf(file)
