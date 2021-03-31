library(HessianScreening)

#d <- readRDS(file.path("data", paste0("dorothea", ".rds")))
d <- generateDesign(100, 10)
X <- d$X
y <- d$y
n <- nrow(X)
p <- ncol(X)
family <- "gaussian"
tol_gap <- 1e-4
tol_infeas <- 1e-3


sparsity <- 1 - Matrix::nnzero(X) / length(X)
sparsity * min(n, p) / max(n, p) * 10

fit_hessian <- lassoPath(
    X,
    y,
    family = family,
    screening_type = "edpp",
    verbosity = 1,
    tol_gap = tol_gap,
    tol_infeas = tol_infeas
)

fit_hessian$lookahead

fit_working <- lassoPath(
    X,
    y,
    family = family,
    screening_type = "working",
    verbosity = 1,
    tol_gap = tol_gap,
    tol_infeas = tol_infeas,
    line_search = TRUE

fit_edpp <- lassoPath(
    X,
    y,
    family = family,
    screening_type = "edpp",
    verbosity = 1,
    tol_gap = tol_gap,
    tol_infeas = tol_infeas
)
fit_gapsafe <- lassoPath(
    X,
    y,
    family = family,
    screening_type = "gap_safe",
    verbosity = 1,
    tol_gap = tol_gap,
    tol_infeas = tol_infeas
)

# # cat("***************\n")
# cat("hessian:\n")
# cat("full = ", fit$full_time, "\n")
# cat("cd_time = ", sum(fit$cd_time), "\n")
# cat("kkt_time = ", sum(fit$kkt_time), "\n")
# cat("hess_time = ", sum(fit$hess_time), "\n")
# cat("gradcorr_time", sum(fit$gradcorr_time), "\n")
# cat("***************\n")
# cat("working:\n")
# cat("full = ", fit.w$full_time, "\n")
# cat("cd_time = ", sum(fit.w$cd_time), "\n")
# cat("kkt_time = ", sum(fit.w$kkt_time), "\n")


# plot(fit$passes, type = "l")
# lines(fit.w$passes, col = "red")
