#pragma once

#include <RcppArmadillo.h>

using namespace arma;

bool
checkStoppingConditions(const uword step,
                        const uword n,
                        const uword p,
                        const uword n_lambda,
                        const uword n_active,
                        const double lambda,
                        const double lambda_min,
                        const double dev,
                        const double dev_prev,
                        const double null_dev,
                        const std::string screening_type,
                        const uword verbosity)
{
  if (step == 1) {
    // first step; always continue
    return false;
  }

  if (verbosity >= 1) {
    Rprintf("  checking stopping conditions\n");
  }

  double dev_ratio = 1.0 - dev / null_dev;
  double dev_change = 1.0 - dev / dev_prev;

  if (verbosity >= 1) {
    Rprintf(
      "    dev ratio:  %.3f\n    dev change: %.6f\n", dev_ratio, dev_change);
  }

  if (dev_change <= 1e-5) {
    return true;
  }

  if (dev_ratio >= 0.999 || lambda <= lambda_min) {
    return true;
  }

  if (n <= p && n_active >= n) {
    return true;
  }

  if (screening_type != "hessian_adaptive" && step >= n_lambda) {
    return true;
  }

  return false;
}
