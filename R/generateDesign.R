#' Generate Response and Design Matrix
#'
#' @param n Number of observations
#' @param p Number of predictors
#' @param family Response type
#' @param density Density of design matrix
#' @param rho Correlation between predictors (corresponding to choice in `rho_type`)
#' @param s Number of signals (nonzero coefficients)
#' @param rho_type Type of correlation between predictors. `"constant"` gives uniform
#'   correlation between all predictors. `"auto"` means that predictors are
#'   "autocorrelated" with adjacent predictors by the amount specified in `rho`.
#' @param beta_type Type of coefficients to generate
#' @param snr Signal-to-noise ratio
#'
#' @return A list with `X`: the design matrix, and `y`: the response.
#' @export
generateDesign <- function(n,
                           p,
                           family = c("gaussian", "binomial"),
                           density = 1,
                           rho = 0,
                           s = 5,
                           rho_type = c("constant", "auto"),
                           beta_type = 1,
                           snr = 2) {
  family <- match.arg(family)
  rho_type <- match.arg(rho_type)

  s <- min(s, p)
  beta <- double(p)

  if (density != 1 && rho > 0) {
    stop("when density != 1, 'rho' must be 0")
  }

  if (beta_type == 1) {
    ind <- round(seq(1, p, length.out = s))
    beta[ind] <- 1
  } else if (beta_type == 2) {
    beta[1:s] <- 1
  } else if (beta_type == 3) {
    beta[1:s] <- seq(10, 0.5, length = s)
  } else if (beta_type == 4) {
    beta[1:s] <- 1
    beta[(s + 1):p] <- 0.5^(1:(p - s))
  }

  if (density == 1) {
    X <- matrix(rnorm(n * p), n)
  } else {
    X <- Matrix::rsparsematrix(n, p, density)
  }

  if (rho != 0) {
    if (rho_type == "auto") {
      inds <- 1:p
      Sigma <- rho^abs(outer(inds, inds, "-"))
      U <- chol(Sigma)

      X <- X %*% t(U)
      sigma <- sqrt((t(beta) %*% Sigma %*% beta) / snr)
    } else if (rho_type == "constant") {
      X <- matrix(rnorm(n), n, p) * sqrt(rho) +
        sqrt(1 - rho) * matrix(rnorm(n*p), n, p)
      sigma <- sqrt((rho * sum(beta)^2 + (1 - rho) * sum(beta^2)) / snr)
    }
  } else {
    Sigma <- diag(p)
    sigma <- sqrt((t(beta) %*% Sigma %*% beta) / snr)
  }

  y <- as.vector(X %*% beta) + as.vector(sigma) * rnorm(n)

  if (family == "binomial") {
    y <- (sign(y) + 1) / 2
  }

  list(X = X, y = y)
}
