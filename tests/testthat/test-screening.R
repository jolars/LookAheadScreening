test_that("screening methods work", {
  set.seed(1014)

  n <- 50
  p <- 200

  for (density in c(1, 0.5)) {
    d <- generateDesign(n, p)

    X <- d$X
    y <- d$y

    fit_gap <- lassoPath(X, y, screening_type = "gap_safe")

    for (screening_type in c("gap_safe_lookahead")) {
      fit <- lassoPath(X,
                       y,
                       screening_type = screening_type,
                       verbosity = 0,
                       force_kkt_check = TRUE)

      steps <- 1:min(length(fit$lambda), length(fit_gap$lambda))
      expect_equal(fit_gap$dev[steps], fit$dev[steps],
                    tolerance = 1e-3)

      expect_equal(sum(fit$violations), 0)
    }
  }
})
