// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rdirichletC
Rcpp::NumericVector rdirichletC(int n);
RcppExport SEXP _SAMEclustering_rdirichletC(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichletC(n));
    return rcpp_result_gen;
END_RCPP
}
// tableC
std::map<int, int> tableC(Rcpp::NumericVector x);
RcppExport SEXP _SAMEclustering_tableC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(tableC(x));
    return rcpp_result_gen;
END_RCPP
}
// EM
Rcpp::List EM(const Rcpp::NumericMatrix& data, const int M);
RcppExport SEXP _SAMEclustering_EM(SEXP dataSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EM(data, M));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SAMEclustering_rdirichletC", (DL_FUNC) &_SAMEclustering_rdirichletC, 1},
    {"_SAMEclustering_tableC", (DL_FUNC) &_SAMEclustering_tableC, 1},
    {"_SAMEclustering_EM", (DL_FUNC) &_SAMEclustering_EM, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SAMEclustering(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
