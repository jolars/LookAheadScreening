#pragma once

#include <RcppArmadillo.h>

using namespace arma;

void
rescaleCoefficients(mat& betas,
                    const vec& X_mean,
                    const vec& X_sd,
                    const double y_mean)
{
  const uword p = betas.n_rows;

  for (uword j = 0; j < p; ++j) {
    betas.row(j) /= X_sd(j);
  }
}