#include <RcppArmadillo.h>

#include "setupModel.h"
#include "model.h"
#include "colNormsSquared.h"
#include "gaussian.h"
#include "checkStoppingConditions.h"
#include "kktCheck.h"
#include "rescaleCoefficients.h"
#include "screenPredictors.h"
#include "standardize.h"

using namespace arma;
using namespace Rcpp;

template<typename T>
Rcpp::List
lassoPath(T& X,
          vec& y,
          const bool standardize,
          const std::string screening_type,
          const uword path_length,
          const uword maxit,
          const double tol_infeas,
          const double tol_gap,
          const bool check_kkt,
          const uword verbosity)
{
  const uword n = X.n_rows;
  const uword p = X.n_cols;

  umat lookaheads(X.n_cols, 0);

  vec beta(p, fill::zeros);
  mat betas(p, 0, fill::zeros);
  vec Xbeta(n, fill::zeros);
  vec residual(n, fill::zeros);
  vec c(p);
  vec c_pred(p);
  vec c_grad(p);

  // standardize predictors and response
  vec X_mean = zeros<vec>(p);
  vec X_sd = ones<vec>(p);

  if (standardize) {
    standardizeX(X_mean, X_sd, X);
  }

  const vec X_mean_scaled = X_mean / X_sd;
  const double y_center = mean(y);

  vec X_norms_squared(p);

  if (!standardize) {
    X_norms_squared = colNormsSquared(X);
  } else {
    X_norms_squared.fill(n);
  }

  auto model = setupModel("gaussian",
                          y,
                          beta,
                          residual,
                          Xbeta,
                          c,
                          X_mean_scaled,
                          X_norms_squared,
                          n,
                          p,
                          standardize);

  model->standardizeY();
  model->updateResidual();
  model->updateCorrelation(X, regspace<uvec>(0, p - 1));

  const double lambda_min_ratio = n < p ? 0.01 : 1e-4;
  const double lambda_max = max(abs(c));
  const double lambda_min = lambda_max * lambda_min_ratio;

  const vec lambdas =
    exp(linspace(log(lambda_max), log(lambda_min), path_length));

  double lambda = lambda_max;

  std::vector<double> primals;
  std::vector<double> duals;
  std::vector<double> devs;
  std::vector<double> dev_ratios;

  std::vector<uword> n_active;
  std::vector<uword> n_new_active;
  std::vector<uword> n_passes;
  std::vector<uword> n_refits;
  std::vector<double> n_screened;
  std::vector<uword> n_violations;

  uvec active(p, fill::zeros);

  uvec first_active = find(abs(c) == lambda_max);

  active(first_active(0)) = true;

  uvec active_prev = active;
  uvec ever_active = active;
  uvec screened = active;

  uvec lookahead(p, fill::zeros);
  uvec lookahead_disabled(p, fill::zeros);

  uvec violations(p, fill::zeros);

  uvec active_set = find(active);
  uvec active_set_prev = active_set;

  vec s(p, fill::zeros);
  s(active_set) = sign(c(active_set));

  const double null_dev = model->deviance();
  double dev = null_dev;

  std::vector<double> it_times;
  std::vector<double> cd_times;

  const double null_primal = model->primal(lambda_max, active_set);

  wall_clock timer;
  timer.tic();

  double full_time = timer.toc();

  uword i = 0;

  while (true) {
    lambda = lambdas(i);

    double it_time = timer.toc();

    if (verbosity >= 1) {
      Rprintf("step: %i, lambda: %.2f\n", i + 1, lambda);
    }

    vec beta_prev = beta;
    double dev_prev = dev;
    uword n_passes_i_sum = 0;
    uword n_violations_i = 0;
    uword n_refits_i = 0;
    bool first_run = true;

    double cd_time = 0;

    if (verbosity >= 1) {
      Rprintf("  running coordinate descent\n");
    }

    double t0 = timer.toc();

    // ever-active warm start
    auto [primal_value, dual_value, duality_gap, n_passes_i, avg_screened] =
      model->fit(ever_active,
                 X,
                 X_norms_squared,
                 lambda,
                 lambda_max,
                 null_primal,
                 "working",
                 true,
                 maxit,
                 tol_gap,
                 tol_infeas,
                 verbosity);

    std::tie(primal_value, dual_value, duality_gap, n_passes_i, avg_screened) =
      model->fit(screened,
                 X,
                 X_norms_squared,
                 lambda,
                 lambda_max,
                 null_primal,
                 screening_type,
                 false,
                 maxit,
                 tol_gap,
                 tol_infeas,
                 verbosity);

    cd_time += timer.toc() - t0;

    n_passes_i_sum += n_passes_i;
    n_screened.push_back(avg_screened);

    t0 = timer.toc();

    if (check_kkt) {
      violations.fill(false);
      const uvec check_set = find(screened == false);
      model->updateCorrelation(X, check_set);
      kktCheck(violations, screened, c, check_set, lambda);
    }

    n_violations_i += sum(violations);

    if (any(violations)) {
      Rcpp::stop("violations ocurred!");
    }

    duals.emplace_back(dual_value);
    primals.emplace_back(primal_value);
    n_passes.emplace_back(n_passes_i_sum);
    n_refits.emplace_back(n_refits_i);
    n_violations.emplace_back(n_violations_i);
    cd_times.emplace_back(cd_time);

    if (i > 0) {
      active = beta != 0;
      active_set = find(active);
      s.zeros();
      s(active_set) = sign(c(active_set));
    }

    dev = model->deviance();
    devs.emplace_back(dev);
    dev_ratios.emplace_back(1.0 - dev / null_dev);

    uword new_active = setDiff(active_set, active_set_prev).n_elem;
    ever_active(active_set).fill(true);
    n_active.emplace_back(active_set.n_elem);
    n_new_active.emplace_back(new_active);

    betas.insert_cols(betas.n_cols, beta);

    if (verbosity >= 1) {
      Rprintf("  active: %i, new active: %i\n", active_set.n_elem, new_active);
    }

    bool stop_path = checkStoppingConditions(i + 1,
                                             n,
                                             p,
                                             path_length,
                                             active_set.n_elem,
                                             lambda,
                                             lambda_min,
                                             dev,
                                             dev_prev,
                                             null_dev,
                                             screening_type,
                                             verbosity);

    if (stop_path) {
      it_times.emplace_back(timer.toc() - it_time);
      break;
    }

    double lambda_next = lambdas(i + 1);

    screenPredictors(lookahead,
                     lookahead_disabled,
                     screened,
                     model,
                     screening_type,
                     ever_active,
                     residual,
                     c,
                     c_grad,
                     X,
                     X_norms_squared,
                     X_mean_scaled,
                     y,
                     lambdas,
                     lambda,
                     lambda_next,
                     i,
                     standardize);

    lookaheads.insert_cols(lookaheads.n_cols, lookahead);

    active_set_prev = active_set;
    lambda = lambda_next;

    it_times.emplace_back(timer.toc() - it_time);

    Rcpp::checkUserInterrupt();

    ++i;
  }

  rescaleCoefficients(betas, X_mean, X_sd, y_center);

  full_time = timer.toc() - full_time;

  return List::create(Named("beta") = wrap(betas),
                      Named("lambda") = wrap(lambdas),
                      Named("primals") = wrap(primals),
                      Named("duals") = wrap(duals),
                      Named("dev_ratio") = wrap(dev_ratios),
                      Named("dev") = wrap(devs),
                      Named("violations") = wrap(n_violations),
                      Named("refits") = wrap(n_refits),
                      Named("active") = wrap(n_active),
                      Named("screened") = wrap(n_screened),
                      Named("new_active") = wrap(n_new_active),
                      Named("passes") = wrap(n_passes),
                      Named("lookahead") = wrap(lookaheads),
                      Named("full_time") = wrap(full_time),
                      Named("it_time") = wrap(it_times),
                      Named("cd_time") = wrap(cd_times));
}

// [[Rcpp::export]]
Rcpp::List
lassoPathDense(arma::mat X,
               arma::vec y,
               const bool standardize,
               const std::string screening_type,
               const arma::uword path_length,
               const arma::uword maxit,
               const double tol_infeas,
               const double tol_gap,
               const bool check_kkt,
               const arma::uword verbosity)
{
  return lassoPath(X,
                   y,
                   standardize,
                   screening_type,
                   path_length,
                   maxit,
                   tol_infeas,
                   tol_gap,
                   check_kkt,
                   verbosity);
}

// [[Rcpp::export]]
Rcpp::List
lassoPathSparse(arma::sp_mat X,
                arma::vec y,
                const bool standardize,
                const std::string screening_type,
                const arma::uword path_length,
                const arma::uword maxit,
                const double tol_infeas,
                const double tol_gap,
                const bool check_kkt,
                const arma::uword verbosity)
{
  return lassoPath(X,
                   y,
                   standardize,
                   screening_type,
                   path_length,
                   maxit,
                   tol_infeas,
                   tol_gap,
                   check_kkt,
                   verbosity);
}
