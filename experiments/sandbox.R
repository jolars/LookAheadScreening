library(LookAheadScreening)

d <- readRDS(file.path("data", paste0("e2006-tfidf-test", ".rds")))
#d <- generateDesign(100, 10000, s = 10)
X <- d$X
y <- d$y
n <- nrow(X)
p <- ncol(X)
family <- "gaussian"
tol_gap <- 1e-5
tol_infeas <- 1e-4

fit <- lassoPath(
    X,
    y,
    screening_type = "gap_safe",
    standardize = TRUE,
    verbosity = 1,
    tol_gap = tol_gap,
    tol_infeas = tol_infeas,
    path_length = 100
)

fit_lookahead <- lassoPath(
    X,
    y,
    screening_type = "gap_safe_lookahead",
    force_kkt_check = FALSE,
    standardize = TRUE,
    verbosity = 1,
    tol_gap = tol_gap,
    tol_infeas = tol_infeas,
    path_length = 100
)

fit$full_time
fit_lookahead$full_time

fit_lookahead$lookahead
fit_lookahead$violations
fit$violations

any(fit_lookahead$violations > 0)


fit_lookahead$lookahead[1:30, 1:39]
