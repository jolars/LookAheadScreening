// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// lassoPathDense
Rcpp::List lassoPathDense(arma::mat X, arma::vec y, const bool standardize, const std::string screening_type, const arma::uword path_length, const arma::uword maxit, const double tol_infeas, const double tol_gap, const bool check_kkt, const arma::uword verbosity);
RcppExport SEXP _LookAheadScreening_lassoPathDense(SEXP XSEXP, SEXP ySEXP, SEXP standardizeSEXP, SEXP screening_typeSEXP, SEXP path_lengthSEXP, SEXP maxitSEXP, SEXP tol_infeasSEXP, SEXP tol_gapSEXP, SEXP check_kktSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool >::type standardize(standardizeSEXP);
    Rcpp::traits::input_parameter< const std::string >::type screening_type(screening_typeSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type path_length(path_lengthSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type tol_infeas(tol_infeasSEXP);
    Rcpp::traits::input_parameter< const double >::type tol_gap(tol_gapSEXP);
    Rcpp::traits::input_parameter< const bool >::type check_kkt(check_kktSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(lassoPathDense(X, y, standardize, screening_type, path_length, maxit, tol_infeas, tol_gap, check_kkt, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// lassoPathSparse
Rcpp::List lassoPathSparse(arma::sp_mat X, arma::vec y, const bool standardize, const std::string screening_type, const arma::uword path_length, const arma::uword maxit, const double tol_infeas, const double tol_gap, const bool check_kkt, const arma::uword verbosity);
RcppExport SEXP _LookAheadScreening_lassoPathSparse(SEXP XSEXP, SEXP ySEXP, SEXP standardizeSEXP, SEXP screening_typeSEXP, SEXP path_lengthSEXP, SEXP maxitSEXP, SEXP tol_infeasSEXP, SEXP tol_gapSEXP, SEXP check_kktSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool >::type standardize(standardizeSEXP);
    Rcpp::traits::input_parameter< const std::string >::type screening_type(screening_typeSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type path_length(path_lengthSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type tol_infeas(tol_infeasSEXP);
    Rcpp::traits::input_parameter< const double >::type tol_gap(tol_gapSEXP);
    Rcpp::traits::input_parameter< const bool >::type check_kkt(check_kktSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(lassoPathSparse(X, y, standardize, screening_type, path_length, maxit, tol_infeas, tol_gap, check_kkt, verbosity));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LookAheadScreening_lassoPathDense", (DL_FUNC) &_LookAheadScreening_lassoPathDense, 10},
    {"_LookAheadScreening_lassoPathSparse", (DL_FUNC) &_LookAheadScreening_lassoPathSparse, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_LookAheadScreening(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
