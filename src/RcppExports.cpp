// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// log_sum_exp
double log_sum_exp(double u, double v);
RcppExport SEXP _juniper0_log_sum_exp(SEXP uSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum_exp(u, v));
    return rcpp_result_gen;
END_RCPP
}
// log_subtract_exp
double log_subtract_exp(double u, double v);
RcppExport SEXP _juniper0_log_subtract_exp(SEXP uSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(log_subtract_exp(u, v));
    return rcpp_result_gen;
END_RCPP
}
// log_sum_exp_vec
double log_sum_exp_vec(NumericVector w);
RcppExport SEXP _juniper0_log_sum_exp_vec(SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(log_sum_exp_vec(w));
    return rcpp_result_gen;
END_RCPP
}
// wbar
NumericVector wbar(double tinf, double dateT, double rOff, double pOff, double pi, double shGen, double scGen, double shSam, double scSam, double delta_t);
RcppExport SEXP _juniper0_wbar(SEXP tinfSEXP, SEXP dateTSEXP, SEXP rOffSEXP, SEXP pOffSEXP, SEXP piSEXP, SEXP shGenSEXP, SEXP scGenSEXP, SEXP shSamSEXP, SEXP scSamSEXP, SEXP delta_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tinf(tinfSEXP);
    Rcpp::traits::input_parameter< double >::type dateT(dateTSEXP);
    Rcpp::traits::input_parameter< double >::type rOff(rOffSEXP);
    Rcpp::traits::input_parameter< double >::type pOff(pOffSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type shGen(shGenSEXP);
    Rcpp::traits::input_parameter< double >::type scGen(scGenSEXP);
    Rcpp::traits::input_parameter< double >::type shSam(shSamSEXP);
    Rcpp::traits::input_parameter< double >::type scSam(scSamSEXP);
    Rcpp::traits::input_parameter< double >::type delta_t(delta_tSEXP);
    rcpp_result_gen = Rcpp::wrap(wbar(tinf, dateT, rOff, pOff, pi, shGen, scGen, shSam, scSam, delta_t));
    return rcpp_result_gen;
END_RCPP
}
// probTTree
double probTTree(NumericMatrix ttree, double rOff, double pOff, double pi, double shGen, double scGen, double shSam, double scSam, double dateT, NumericVector wbar0, double delta_t);
RcppExport SEXP _juniper0_probTTree(SEXP ttreeSEXP, SEXP rOffSEXP, SEXP pOffSEXP, SEXP piSEXP, SEXP shGenSEXP, SEXP scGenSEXP, SEXP shSamSEXP, SEXP scSamSEXP, SEXP dateTSEXP, SEXP wbar0SEXP, SEXP delta_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ttree(ttreeSEXP);
    Rcpp::traits::input_parameter< double >::type rOff(rOffSEXP);
    Rcpp::traits::input_parameter< double >::type pOff(pOffSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type shGen(shGenSEXP);
    Rcpp::traits::input_parameter< double >::type scGen(scGenSEXP);
    Rcpp::traits::input_parameter< double >::type shSam(shSamSEXP);
    Rcpp::traits::input_parameter< double >::type scSam(scSamSEXP);
    Rcpp::traits::input_parameter< double >::type dateT(dateTSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type wbar0(wbar0SEXP);
    Rcpp::traits::input_parameter< double >::type delta_t(delta_tSEXP);
    rcpp_result_gen = Rcpp::wrap(probTTree(ttree, rOff, pOff, pi, shGen, scGen, shSam, scSam, dateT, wbar0, delta_t));
    return rcpp_result_gen;
END_RCPP
}
// dprop
std::vector<double> dprop(std::vector<double> x, double mu, bool LOG);
RcppExport SEXP _juniper0_dprop(SEXP xSEXP, SEXP muSEXP, SEXP LOGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< bool >::type LOG(LOGSEXP);
    rcpp_result_gen = Rcpp::wrap(dprop(x, mu, LOG));
    return rcpp_result_gen;
END_RCPP
}
// dprop_bounded
std::vector<double> dprop_bounded(std::vector<double> x, std::vector<double> N, double mu, bool LOG);
RcppExport SEXP _juniper0_dprop_bounded(SEXP xSEXP, SEXP NSEXP, SEXP muSEXP, SEXP LOGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< bool >::type LOG(LOGSEXP);
    rcpp_result_gen = Rcpp::wrap(dprop_bounded(x, N, mu, LOG));
    return rcpp_result_gen;
END_RCPP
}
// ppropBoundedHelper
double ppropBoundedHelper(double x, double N, double mu);
RcppExport SEXP _juniper0_ppropBoundedHelper(SEXP xSEXP, SEXP NSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(ppropBoundedHelper(x, N, mu));
    return rcpp_result_gen;
END_RCPP
}
// pprop_bounded
std::vector<double> pprop_bounded(double x, std::vector<double> N, double mu);
RcppExport SEXP _juniper0_pprop_bounded(SEXP xSEXP, SEXP NSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(pprop_bounded(x, N, mu));
    return rcpp_result_gen;
END_RCPP
}
// pprop
std::vector<double> pprop(std::vector<double> x, double mu, bool LOG);
RcppExport SEXP _juniper0_pprop(SEXP xSEXP, SEXP muSEXP, SEXP LOGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< bool >::type LOG(LOGSEXP);
    rcpp_result_gen = Rcpp::wrap(pprop(x, mu, LOG));
    return rcpp_result_gen;
END_RCPP
}
// sprop
std::vector<double> sprop(std::vector<double> x, double mu, bool LOG);
RcppExport SEXP _juniper0_sprop(SEXP xSEXP, SEXP muSEXP, SEXP LOGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< bool >::type LOG(LOGSEXP);
    rcpp_result_gen = Rcpp::wrap(sprop(x, mu, LOG));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_juniper0_log_sum_exp", (DL_FUNC) &_juniper0_log_sum_exp, 2},
    {"_juniper0_log_subtract_exp", (DL_FUNC) &_juniper0_log_subtract_exp, 2},
    {"_juniper0_log_sum_exp_vec", (DL_FUNC) &_juniper0_log_sum_exp_vec, 1},
    {"_juniper0_wbar", (DL_FUNC) &_juniper0_wbar, 10},
    {"_juniper0_probTTree", (DL_FUNC) &_juniper0_probTTree, 11},
    {"_juniper0_dprop", (DL_FUNC) &_juniper0_dprop, 3},
    {"_juniper0_dprop_bounded", (DL_FUNC) &_juniper0_dprop_bounded, 4},
    {"_juniper0_ppropBoundedHelper", (DL_FUNC) &_juniper0_ppropBoundedHelper, 3},
    {"_juniper0_pprop_bounded", (DL_FUNC) &_juniper0_pprop_bounded, 3},
    {"_juniper0_pprop", (DL_FUNC) &_juniper0_pprop, 3},
    {"_juniper0_sprop", (DL_FUNC) &_juniper0_sprop, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_juniper0(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
