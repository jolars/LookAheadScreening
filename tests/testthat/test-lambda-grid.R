library(HessianScreening)

test_that("lambda grid calculations are correct", {
  for (p in c(10, 100)) {
    set.seed(1110)

    n <- 50

    X <- matrix(rnorm(n * p), n, p)
    beta <- rnorm(p)
    y <- X %*% beta + rnorm(n)

    X <- scale(X)
    y <- y - mean(y)

    family <- "gaussian"

    res_work <- lassoPath(X, y, family = family)
    res_glmn <- glmnet::glmnet(X, y, intercept = FALSE)

    n_lambda <- min(length(res_work$lambda), length(res_glmn$lambda))

    expect_equal(res_work$lambda[1:n_lambda] / n, res_glmn$lambda[1:n_lambda])
  }
})
