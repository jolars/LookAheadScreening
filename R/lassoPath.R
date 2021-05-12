#' Lasso Path with Hessian Screening Rules
#'
#' @param X The predictor matrix
#' @param y The reponse vector
#' @param standardize Whether to standardize the predictors
#' @param screening_type Screening rule
#' @param path_length The (desired) length of the lasso path
#' @param maxit Maximum number of iterations for Coordinate Descent loop
#' @param tol_infeas Tolerance threshold for maximum infeasibility
#' @param tol_gap Tolerance threshold for duality gap
#' @param check_kkt Whether to force KKT checks
#' @param verbosity Controls the level of verbosity. 0 = no output.
#'
#' @export
lassoPath <- function(X,
                      y,
                      standardize = TRUE,
                      screening_type = c(
                        "gap_safe",
                        "gap_safe_lookahead"
                      ),
                      path_length = 100L,
                      maxit = 1e5,
                      tol_infeas = 1e-4,
                      tol_gap = 1e-5,
                      check_kkt = FALSE,
                      verbosity = 0) {
  screening_type <- match.arg(screening_type)

  sparse <- inherits(X, "sparseMatrix")

  n <- nrow(X)
  p <- ncol(X)

  if (sparse) {
    X <- methods::as(X, "dgCMatrix")

    lassoPathSparse(
      X,
      y,
      standardize,
      screening_type,
      path_length,
      maxit,
      tol_infeas,
      tol_gap,
      check_kkt,
      verbosity
    )
  } else {
    lassoPathDense(
      X,
      y,
      standardize,
      screening_type,
      path_length,
      maxit,
      tol_infeas,
      tol_gap,
      check_kkt,
      verbosity
    )
  }
}
