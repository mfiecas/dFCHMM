// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// colSums_cpp
arma::vec colSums_cpp(const arma::mat& A);
RcppExport SEXP _dFCHMM_colSums_cpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(colSums_cpp(A));
    return rcpp_result_gen;
END_RCPP
}
// rowSums_cpp
arma::vec rowSums_cpp(const arma::mat& A);
RcppExport SEXP _dFCHMM_rowSums_cpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(rowSums_cpp(A));
    return rcpp_result_gen;
END_RCPP
}
// logvec_c
arma::vec logvec_c(const arma::vec& x);
RcppExport SEXP _dFCHMM_logvec_c(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logvec_c(x));
    return rcpp_result_gen;
END_RCPP
}
// factorial_cpp
unsigned long factorial_cpp(unsigned n);
RcppExport SEXP _dFCHMM_factorial_cpp(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(factorial_cpp(n));
    return rcpp_result_gen;
END_RCPP
}
// covmat_c
arma::mat covmat_c(arma::mat Xt, arma::vec mu, arma::vec u);
RcppExport SEXP _dFCHMM_covmat_c(SEXP XtSEXP, SEXP muSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xt(XtSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(covmat_c(Xt, mu, u));
    return rcpp_result_gen;
END_RCPP
}
// standardize_rows_cpp
arma::mat standardize_rows_cpp(arma::mat A);
RcppExport SEXP _dFCHMM_standardize_rows_cpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(standardize_rows_cpp(A));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrm_arma
arma::vec dmvnrm_arma(const arma::mat& x, const arma::rowvec& mean, const arma::mat& sigma, bool logd);
RcppExport SEXP _dFCHMM_dmvnrm_arma(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm_arma(x, mean, sigma, logd));
    return rcpp_result_gen;
END_RCPP
}
// eye_cpp
arma::mat eye_cpp(int P);
RcppExport SEXP _dFCHMM_eye_cpp(SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(eye_cpp(P));
    return rcpp_result_gen;
END_RCPP
}
// logsumexp_cpp
double logsumexp_cpp(arma::vec X);
RcppExport SEXP _dFCHMM_logsumexp_cpp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(logsumexp_cpp(X));
    return rcpp_result_gen;
END_RCPP
}
// forward_backward_cpp
List forward_backward_cpp(const arma::mat& Xt, const arma::mat& B, const arma::mat& mu, const arma::cube& Sigma, const arma::vec& delta);
RcppExport SEXP _dFCHMM_forward_backward_cpp(SEXP XtSEXP, SEXP BSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Xt(XtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(forward_backward_cpp(Xt, B, mu, Sigma, delta));
    return rcpp_result_gen;
END_RCPP
}
// hmm_cpp
List hmm_cpp(const arma::cube& Xt, int num_states, const arma::cube& Sigma_init, int maxiter, double tol, bool verbose);
RcppExport SEXP _dFCHMM_hmm_cpp(SEXP XtSEXP, SEXP num_statesSEXP, SEXP Sigma_initSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type Xt(XtSEXP);
    Rcpp::traits::input_parameter< int >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Sigma_init(Sigma_initSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(hmm_cpp(Xt, num_states, Sigma_init, maxiter, tol, verbose));
    return rcpp_result_gen;
END_RCPP
}
