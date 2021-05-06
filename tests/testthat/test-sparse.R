test_that("sparse and dense methods are equivalent", {
  n <- 100
  p <- 5

  d <- generateDesign(n, p, density = 0.5)
  X <- d$X
  y <- d$y
  for (standardize in c(TRUE, FALSE)) {
    fit_dense <- lassoPath(as.matrix(X), y, standardize = standardize)
    fit_sparse <- lassoPath(X, y, standardize = standardize)

    expect_equal(fit_dense$lambda, fit_sparse$lambda)
    expect_equal(fit_dense$dev_ratio, fit_sparse$dev_ratio)
    expect_equal(fit_dense$beta, fit_sparse$beta)
  }
})
