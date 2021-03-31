#pragma once

#include <RcppArmadillo.h>

#include "prox.h"

using namespace arma;

class Model
{
public:
  const std::string family;

  vec& y;
  vec& beta;
  vec& residual;
  vec& Xbeta;
  vec& c;

  const vec& X_mean_scaled;
  const vec& X_norms_squared;

  const uword n;
  const uword p;
  const bool standardize;

  double dual_scale{ 0 };

  Model(const std::string family,
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
    : family(family)
    , y(y)
    , beta(beta)
    , residual(residual)
    , Xbeta(Xbeta)
    , c(c)
    , X_mean_scaled(X_mean_scaled)
    , X_norms_squared(X_norms_squared)
    , n(n)
    , p(p)
    , standardize(standardize)
  {}

  virtual ~Model() = default;

  virtual double primal(const double lambda, const uvec& screened_set) = 0;

  virtual double dual() = 0;

  virtual double scaledDual(const double lambda) = 0;

  virtual double deviance() = 0;

  virtual double hessianTerm(const mat& X, const uword j) = 0;

  virtual double hessianTerm(const sp_mat& X, const uword j) = 0;

  virtual void updateResidual() = 0;

  virtual void adjustResidual(const mat& X,
                              const uword j,
                              const double beta_diff) = 0;

  virtual void adjustResidual(const sp_mat& X,
                              const uword j,
                              const double beta_diff) = 0;

  void updateLinearPredictor(const mat& X, const uvec& ind)
  {
    Xbeta = X.cols(ind) * beta(ind);
  }

  void updateLinearPredictor(const sp_mat& X, const uvec& ind)
  {
    Xbeta = X.cols(ind) * beta(ind);

    if (standardize)
      Xbeta -= dot(beta(ind), X_mean_scaled(ind));
  }

  void updateCorrelation(const mat& X, const uvec& ind)
  {
    for (auto&& j : ind) {
      c(j) = dot(X.unsafe_col(j), residual);
    }
  }

  void updateCorrelation(const sp_mat& X, const uvec& ind)
  {
    for (auto&& j : ind) {
      c(j) = dot(X.col(j), residual);
    }

    if (standardize) {
      c(ind) -= X_mean_scaled(ind) * sum(residual);
    }
  }

  inline void updateCorrelation(const mat& X, const uword& j)
  {
    c(j) = dot(X.unsafe_col(j), residual);
  }

  inline void updateCorrelation(const sp_mat& X, const uword& j)
  {
    c(j) = dot(X.col(j), residual);

    if (standardize) {
      c(j) -= X_mean_scaled(j) * accu(residual);
    }
  }

  virtual void standardizeY() = 0;

  virtual double safeScreeningRadius(const double duality_gap,
                                     const double lambda) = 0;

  template<typename T>
  void safeScreening(uvec& screened,
                     uvec& screened_set,
                     const T& X,
                     const vec& XTcenter,
                     const double r_screen)
  {
    for (auto&& j : screened_set) {
      double r_normX_j = r_screen * std::sqrt(X_norms_squared(j));

      if (r_normX_j >= 1) {
        continue;
      }

      if (std::abs(XTcenter(j)) + r_normX_j + std::sqrt(datum::eps) < 1) {
        // predictor must be zero; update residual and remove from screened set
        if (beta(j) != 0) {
          adjustResidual(X, j, -beta(j));
          beta(j) = 0;
        }

        screened(j) = false;
      }
    }

    screened_set = find(screened);
  }

  template<typename T>
  std::tuple<double, double, double, uword, double> fit(
    uvec& screened,
    const T& X,
    const vec& X_norms_squared,
    const double lambda,
    const double lambda_max,
    const double null_primal,
    const std::string screening_type,
    const bool first_run,
    const uword step,
    const uword maxit,
    const double tol_gap,
    const double tol_infeas,
    const bool line_search,
    const uword verbosity)
  {
    const uword p = X.n_cols;

    if (screening_type == "gap_safe" && !first_run)
      screened.fill(true);

    uvec screened_set = find(screened);

    double primal_value = primal(lambda, screened_set);
    double dual_value = dual();
    double duality_gap = primal_value - dual_value;

    // line search parameters
    const double a = 0.1;
    const double b = 0.5;

    vec XTcenter(p);

    vec t(p, fill::ones); // learning rates

    const uword screen_interval = 10;

    double n_screened = 0;
    uword it = 0;

    if (!screened_set.is_empty()) {
      updateLinearPredictor(X, screened_set);
      updateResidual();

      while (it < maxit) {
        it++;

        if (verbosity >= 2) {
          Rprintf("    iter: %i\n", it);
        }

        if (screening_type == "gap_safe" && (it - 1) % screen_interval == 0) {
          if (it > 1) {
            updateLinearPredictor(X, screened_set);
            updateResidual();
          }
          updateCorrelation(X, screened_set);

          dual_scale = std::max(lambda, max(abs(c)));

          primal_value = primal(lambda, screened_set);
          dual_value = scaledDual(lambda);
          duality_gap = primal_value - dual_value;

          double r_screen{ 0 };

          if (screening_type == "gap_safe") {
            XTcenter = c / dual_scale;
            r_screen = safeScreeningRadius(std::max(duality_gap, 0.0), lambda);
          }

          safeScreening(screened, screened_set, X, XTcenter, r_screen);
        }

        n_screened += screened_set.n_elem;

        for (auto&& j : screened_set) {
          updateCorrelation(X, j);
          double hess_j = hessianTerm(X, j);

          double beta_j_old = beta(j);
          double v =
            prox(beta_j_old + c(j) / hess_j, lambda / hess_j) - beta(j);

          if (v != 0) {
            if (family == "binomial" && line_search) {
              // line search (see J. D. Lee, Y. Sun, and M. A. Saunders,
              // “Proximal Newton-type methods for minimizing composite
              // functions,” arXiv:1206.1623 [cs, math, stat], Mar. 2014,
              // Accessed: Jan. 12, 2020. [Online]. Available:
              // http://arxiv.org/abs/1206.1623)

              double primal_value_old = primal(lambda, screened_set);

              while (true) {
                double beta_j_prev = beta(j);

                beta(j) = beta_j_old + t(j) * v;
                adjustResidual(X, j, beta(j) - beta_j_prev);

                primal_value = primal(lambda, screened_set);

                double eta = -c(j) * v + lambda * (std::abs(beta_j_old + v) -
                                                   std::abs(beta_j_old));

                if (primal_value * (1 - std::sqrt(datum::eps)) <=
                    primal_value_old + a * t(j) * eta) {
                  break;
                } else {
                  t(j) *= b;
                }

                Rcpp::checkUserInterrupt();
              }
            } else {
              beta(j) = beta_j_old + v;
              adjustResidual(X, j, beta(j) - beta_j_old);
            }
          }
        }

        primal_value = primal(lambda, screened_set);
        dual_value = dual();

        duality_gap = primal_value - dual_value;

        if (verbosity >= 2)
          Rprintf("      primal: %f, dual: %f, duality gap: %f\n",
                  primal_value,
                  dual_value,
                  duality_gap / null_primal);

        if (std::abs(duality_gap) <= tol_gap * null_primal) {
          updateCorrelation(X, screened_set);

          double infeas = lambda > 0 ? max(abs(c(screened_set)) - lambda) : 0;

          if (verbosity >= 2)
            Rprintf("      infeasibility: %f\n", infeas / lambda_max);

          if (infeas <= lambda_max * tol_infeas)
            break;
        }

        Rcpp::checkUserInterrupt();
      }
    } else {
      beta.zeros();
    }

    double avg_screened = n_screened / it;

    return { primal_value, dual_value, duality_gap, it, avg_screened };
  }
};
