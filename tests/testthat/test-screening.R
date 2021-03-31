test_that("screening methods work", {
  set.seed(1014)

  n <- 50
  p <- 200

  for (density in c(1, 0.5)) {

    for (family in c("gaussian", "binomial")) {
      d <- generateDesign(n, p, family = family)

      X <- d$X
      y <- d$y

      fit_work <- lassoPath(X, y, family = family, screening_type = "working")

      for (screening_type in c("hessian", "working", "strong", "gap_safe")) {
        if (family == "gaussian" && screening_type == "edpp") {
          next
        }

        fit <- lassoPath(X,
                         y,
                         family = family,
                         screening_type = screening_type,
                         verbosity = 0,
                         force_kkt_check = TRUE)

        steps <- 1:min(length(fit$lambda), length(fit_work$lambda))
        expect_equal(fit_work$dev[steps], fit$dev[steps],
                     tolerance = 1e-3)

        if (screening_type %in% c("gap_safe")) {
          expect_equal(sum(fit$violations), 0)
        }
      }
    }
  }
})
