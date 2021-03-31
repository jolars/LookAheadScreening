test_that("compare against glmnet with real data", {
  suppressMessages(library(glmnet))

  datalist <- c("heart", "leukemia")

  for (dataset in datalist) {
    data(list = list(dataset))
    d <- get(dataset)
    X <- d$X
    y <- d$y

    fit_ours <- lassoPath(X,
      y,
      "binomial",
      screening_type = "working",
      standardize = FALSE
    )

    lambda <- fit_ours$lambda / nrow(X)

    fit_gnet <- glmnet(X,
      y,
      "binomial",
      lambda = lambda,
      intercept = FALSE,
      standardize = FALSE
    )

    ours_dev <- fit_ours$dev_ratio
    gnet_dev <- fit_gnet$dev.ratio

    dev_diff <- pmax(gnet_dev - ours_dev, 0)

    expect_equal(dev_diff, rep(0, length(dev_diff)), tolerance = 1e-3)
  }
})
