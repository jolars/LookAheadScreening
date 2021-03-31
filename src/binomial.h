#pragma once

#include <RcppArmadillo.h>

#include "prox.h"
#include "utils.h"

using namespace arma;

class Binomial : public Model
{
public:
  vec expXbeta;
  vec pr;
  vec w;

  const double p_min = 1e-5;
  const double p_max = 1 - p_min;

  Binomial(const std::string family,
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
    , expXbeta(y.n_elem, fill::zeros)
    , pr(y.n_elem, fill::zeros)
    , w(y.n_elem, fill::zeros)
  {}

  double primal(const double lambda, const uvec& screened_set)
  {
    return -sum(y % Xbeta - log1p(expXbeta)) +
           lambda * norm(beta(screened_set), 1);
  }

  double dual() { return -sum(pr % log(pr) + (1 - pr) % log(1 - pr)); }

  double scaledDual(const double lambda)
  {
    if (dual_scale == 0) {
      return 0;
    } else {
      double alpha = lambda / dual_scale;

      vec prx = clamp(y - alpha * residual, p_min, p_max);

      return -sum(prx % log(prx) + (1 - prx) % log(1 - prx));
    }
  }

  double deviance() { return -2 * sum(y % Xbeta - log1p(expXbeta)); }

  double hessianTerm(const mat& X, const uword j)
  {
    return std::max(dot(square(X.col(j)), w), std::sqrt(datum::eps));
  }

  double hessianTerm(const sp_mat& X, const uword j)
  {
    double out = dot(square(X.col(j)), w);

    if (standardize) {
      out += std::pow(X_mean_scaled(j), 2) * sum(w) -
             2 * dot(X.col(j), w) * X_mean_scaled(j);
    }

    return std::max(out, std::sqrt(datum::eps));
  }

  void updateResidual()
  {
    expXbeta = exp(Xbeta);
    pr = clamp(expXbeta / (1 + expXbeta), p_min, p_max);
    w = pr % (1 - pr);
    residual = y - pr;
  }

  void adjustResidual(const mat& X, const uword j, const double beta_diff)
  {
    Xbeta += X.col(j) * beta_diff;
    updateResidual();
  }

  void adjustResidual(const sp_mat& X, const uword j, const double beta_diff)
  {
    Xbeta += X.col(j) * beta_diff;

    if (standardize)
      Xbeta -= X_mean_scaled(j) * beta_diff;

    updateResidual();
  }

  void standardizeY() {}

  double safeScreeningRadius(const double duality_gap, const double lambda)
  {
    return std::sqrt(2 * duality_gap) / (2 * lambda);
  }
};
;
