// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fasterRndWalk
double fasterRndWalk(IntegerVector gSetIdx, IntegerVector geneRanking, int j, NumericMatrix Ra);
RcppExport SEXP _GSVA_fasterRndWalk(SEXP gSetIdxSEXP, SEXP geneRankingSEXP, SEXP jSEXP, SEXP RaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type gSetIdx(gSetIdxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type geneRanking(geneRankingSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ra(RaSEXP);
    rcpp_result_gen = Rcpp::wrap(fasterRndWalk(gSetIdx, geneRanking, j, Ra));
    return rcpp_result_gen;
END_RCPP
}

RcppExport void ks_matrix_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
RcppExport void matrix_density_R(void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_GSVA_fasterRndWalk", (DL_FUNC) &_GSVA_fasterRndWalk, 4},
    {"ks_matrix_R",      (DL_FUNC) &ks_matrix_R,      10},
    {"matrix_density_R", (DL_FUNC) &matrix_density_R,  7},
    {NULL, NULL, 0}
};

RcppExport void R_init_GSVA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
