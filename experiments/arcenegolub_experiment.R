#rm(list = ls())
#graphics.off()
library(HessianScreening)


#d <- readRDS(file.path("data", paste0("leukemia-train", ".rds")))
#d <- readRDS(file.path("data", paste0("arcene", ".rds")))
#d <- readRDS(file.path("data", paste0("dorothea", ".rds")))
#d <- readRDS(file.path("data", paste0("covtype", ".rds")))
d <- readRDS(file.path("data", paste0("gisette-train", ".rds")))

X <- d$X
y <- d$y
family <- "binomial"
fita <- lassoPath(
    X,
    y,
    family = family,
    screening_type = "hessian",
    hessian_warm_starts = TRUE,
    log_hessian_update_type = "approx",
    gamma = 0.01,
    verify_hessian = FALSE,
    verbosity = 0
)
fitf <- lassoPath(
    X,
    y,
    family = family,
    screening_type = "hessian",
    hessian_warm_starts = TRUE,
    log_hessian_update_type = "full",
    gamma = 0.01,
    verify_hessian = FALSE,
    verbosity = 0
)
fit.w <- lassoPath(
    X,
    y,
    family = family,
    screening_type = "working",
    verbosity = 0
)
cat("***************\n")
cat("hessian:\n")
cat("full = ", fit$full_time, "\n")
cat("cd_time = ", sum(fit$cd_time), "\n")
cat("kkt_time = ", sum(fit$corr_time), "\n")
cat("hess_time = ", sum(fit$hess_time), "\n")
cat("gradcorr_time", sum(fit$gradcorr_time), "\n")
cat("duplicates_time", sum(fit$duplicates_time), "\n")
cat("***************\n")
# cat("working:\n")
cat("full = ", fit.w$full_time, "\n")
cat("cd_time = ", sum(fit.w$cd_time), "\n")
cat("kkt_time = ", sum(fit.w$kkt_time), "\n")
cat("duplicates_time", sum(fit$duplicates_time), "\n")

par(mfrow=c(2,2))
plot(fit$active, type = "l", lty = 2, ylim = c(0, max(c(fit$active, fit$screened))))
lines(fit$screened)

plot(fit$cd_time, type = "l")
lines(fit.w$cd_time, col = "red")

#
plot(fit$passes, type = "l", ylim = c(0, max(c(fit$passes, fit.w$passes))))
lines(fit.w$passes, col = "red")

plot(fit$dev[1:length(fit.w$dev)]-fit.w$dev, type = "l", lwd = 2)

#plot(fit$violations)
