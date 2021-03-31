#pragma once

#include <RcppArmadillo.h>

#include "prox.h"
#include "utils.h"

using namespace arma;

class Gaussian : public Model
{
public:
  Gaussian(const std::string family,
           vec& y,
           vec& beta,
           vec& residual,
           vec& Xbeta,
           vec& c,
           const vec& X_mean_scaled,
           const vec& X_norms_squared,
           const uword n,
           const uword p,
           const bool standardize)
    : Model{ family,          y, beta, residual,   Xbeta, c, X_mean_scaled,
             X_norms_squared, n, p,    standardize }
  {}

  double primal(const double lambda, const uvec& screened_set)
  {
    return 0.5 * std::pow(norm(residual), 2) +
           lambda * norm(beta(screened_set), 1);
  }

  double dual() { return dot(residual, y) - 0.5 * std::pow(norm(residual), 2); }

  double scaledDual(const double lambda)
  {
    if (dual_scale == 0) {
      return 0;
    } else {
      double alpha = lambda / dual_scale;

      return alpha * dot(residual, y) -
             0.5 * std::pow(alpha * norm(residual), 2);
    }
  }

  double deviance() { return std::pow(norm(residual), 2); }

  double hessianTerm(const mat& X, const uword j) { return X_norms_squared(j); }

  double hessianTerm(const sp_mat& X, const uword j)
  {
    return X_norms_squared(j);
  }

  void updateResidual() { residual = y - Xbeta; }

  void adjustResidual(const mat& X, const uword j, const double beta_diff)
  {
    residual -= X.col(j) * beta_diff;
  }

  void adjustResidual(const sp_mat& X, const uword j, const double beta_diff)
  {
    residual -= X.col(j) * beta_diff;

    if (standardize)
      residual += X_mean_scaled(j) * beta_diff;
  }

  void standardizeY() { y -= mean(y); }

  double safeScreeningRadius(const double duality_gap, const double lambda)
  {
    return std::sqrt(duality_gap) / lambda;
  }
};
