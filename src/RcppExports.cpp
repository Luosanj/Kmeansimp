// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// acCtr
NumericVector acCtr(NumericMatrix exdata, NumericMatrix centroids, int i);
RcppExport SEXP _Kmeansimp_acCtr(SEXP exdataSEXP, SEXP centroidsSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type exdata(exdataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centroids(centroidsSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(acCtr(exdata, centroids, i));
    return rcpp_result_gen;
END_RCPP
}
// doClus
NumericVector doClus(NumericMatrix dataset, NumericMatrix centroids);
RcppExport SEXP _Kmeansimp_doClus(SEXP datasetSEXP, SEXP centroidsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dataset(datasetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type centroids(centroidsSEXP);
    rcpp_result_gen = Rcpp::wrap(doClus(dataset, centroids));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Kmeansimp_acCtr", (DL_FUNC) &_Kmeansimp_acCtr, 3},
    {"_Kmeansimp_doClus", (DL_FUNC) &_Kmeansimp_doClus, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_Kmeansimp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}