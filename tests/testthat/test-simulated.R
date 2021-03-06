test_that("gaussian models for simulated data", {
  library(LookAheadScreening)
  library(glmnet)

  grid <- expand.grid(
    np = list(c(100, 5), c(50, 100)),
    density = c(0.2, 1),
    standardize = c(FALSE),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(grid)) {
    set.seed(i)

    g <- grid[i, ]

    np <- g$np[[1]]
    n <- np[1]
    p <- np[2]
    standardize <- g$standardize

    data <- generateDesign(
      n,
      p,
      density = g$density
    )

    X <- data$X
    y <- data$y
    y <- y - mean(y)

    fit_glmnet <- glmnet(
      X,
      y,
      intercept = FALSE,
      standardize = standardize
    )

    fit_ours <- lassoPath(
      X,
      y,
      standardize = standardize
    )

    n_lambda <- min(ncol(fit_glmnet$beta), ncol(fit_ours$beta))

    glmnet_lambda <- fit_glmnet$lambda[1:n_lambda] * n
    ours_lambda <- fit_ours$lambda[1:n_lambda]

    expect_equal(glmnet_lambda, ours_lambda)

    glmnet_beta <- as.matrix(fit_glmnet$beta[, 1:n_lambda])
    ours_beta <- fit_ours$beta[, 1:n_lambda]

    glmnet_devratio <- fit_glmnet$dev.ratio[1:n_lambda]
    ours_devratio <- fit_ours$dev_ratio[1:n_lambda]

    dev_diff <- pmax(glmnet_devratio - ours_devratio, 0)

    expect_equal(dev_diff, rep(0, n_lambda), tolerance = 1e-3)

    if (n > p) {
      expect_equal(ours_beta, glmnet_beta, ignore_attr = TRUE, tolerance = 1e-2)
    }
  }
})
