#pragma once

#include <RcppArmadillo.h>

#include "model.h"

using namespace arma;

template<typename T>
void
screenPredictors(uvec& lookahead,
                 uvec& lookahead_disabled,
                 uvec& screened,
                 const std::unique_ptr<Model>& model,
                 const std::string screening_type,
                 const uvec& ever_active,
                 const vec& residual,
                 const vec& corr,
                 const vec& corr_grad,
                 const T& X,
                 const vec& X_norms_squared,
                 const vec& X_mean_scaled,
                 const vec& y,
                 const vec& lambdas,
                 const double lambda,
                 const double lambda_next,
                 const uword step,
                 const bool standardize)
{
  uword p = X.n_cols;

  if (screening_type == "gap_safe") {
    screened.fill(true);
  } else if (screening_type == "gap_safe_lookahead") {
    // create dual feasible point
    double dual_scale = std::max(lambda, max(abs(corr)));
    vec theta = residual / dual_scale;

    double theta_dot_y = dot(theta, y);
    double theta_dot_theta = dot(theta, theta);
    double beta_norm1 = norm(model->beta, 1);
    double residual_sq_norm2 = std::pow(norm(residual), 2);

    for (uword j = 0; j < p; ++j) {
      if (ever_active(j) || lookahead(j) > step || lookahead_disabled(j))
        continue;

      double a = std::pow(1 - std::abs(corr(j) / dual_scale), 2) -
                 0.5 * theta_dot_theta * X_norms_squared(j);
      double b = X_norms_squared(j) * (theta_dot_y - beta_norm1);
      double c = -0.5 * residual_sq_norm2 * X_norms_squared(j);

      double lambda_star1 = (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a);
      double lambda_star2 = (-b - std::sqrt(b * b - 4 * a * c)) / (2 * a);

      // Rcpp::Rcout << "a = " << a << ", b = " << b << ", c = " << c
      //             << ",lambda_star1: " << lambda_star1
      //             << ", lambda_star2: " << lambda_star2 << std::endl;

      uvec tmp = lambdas > lambda_star1 && lambdas < lambda;

      if (any(tmp)) {
        lookahead(j) = as_scalar(find(tmp, 1, "last"));
      } else {
        // stop considering predictors that are not captured by the rule
        lookahead(j) = 0;
        lookahead_disabled(j) = true;
      }
    }

    // for (uword i = step; i < lambdas.n_elem; ++i) {
    //   if (lambdas(i) < lambda) {
    //     screened.ones();
    //     screened_set = find(screened);

    //     double primal_value = model->primal(lambdas(i), screened_set);
    //     double dual_value = model->scaledDual(lambdas(i));
    //     double duality_gap = primal_value - dual_value;

    //     vec XTcenter = corr / dual_scale;
    //     double r_screen =
    //       model->safeScreeningRadius(std::max(duality_gap, 0.0), lambdas(i));

    //     model->safeScreening(screened, screened_set, X, XTcenter, r_screen);

    //     for (uword j = 0; j < X.n_cols; ++j) {
    //       if (screened(j)) {
    //         lookahead(j) = std::min(lookahead(j), i);
    //       }
    //     }
    //   }
    // }

    screened = lookahead <= step || ever_active;
  }
}
