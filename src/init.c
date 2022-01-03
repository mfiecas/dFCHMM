#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _dFCHMM_colSums_cpp(SEXP);
extern SEXP _dFCHMM_covmat_c(SEXP, SEXP, SEXP);
extern SEXP _dFCHMM_dmvnrm_arma(SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCHMM_eye_cpp(SEXP);
extern SEXP _dFCHMM_factorial_cpp(SEXP);
extern SEXP _dFCHMM_forward_backward_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCHMM_hmm_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCHMM_logsumexp_cpp(SEXP);
extern SEXP _dFCHMM_logvec_c(SEXP);
extern SEXP _dFCHMM_rowSums_cpp(SEXP);
extern SEXP _dFCHMM_standardize_rows_cpp(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_dFCHMM_colSums_cpp",          (DL_FUNC) &_dFCHMM_colSums_cpp,          1},
    {"_dFCHMM_covmat_c",             (DL_FUNC) &_dFCHMM_covmat_c,             3},
    {"_dFCHMM_dmvnrm_arma",          (DL_FUNC) &_dFCHMM_dmvnrm_arma,          4},
    {"_dFCHMM_eye_cpp",              (DL_FUNC) &_dFCHMM_eye_cpp,              1},
    {"_dFCHMM_factorial_cpp",        (DL_FUNC) &_dFCHMM_factorial_cpp,        1},
    {"_dFCHMM_forward_backward_cpp", (DL_FUNC) &_dFCHMM_forward_backward_cpp, 5},
    {"_dFCHMM_hmm_cpp",              (DL_FUNC) &_dFCHMM_hmm_cpp,              6},
    {"_dFCHMM_logsumexp_cpp",        (DL_FUNC) &_dFCHMM_logsumexp_cpp,        1},
    {"_dFCHMM_logvec_c",             (DL_FUNC) &_dFCHMM_logvec_c,             1},
    {"_dFCHMM_rowSums_cpp",          (DL_FUNC) &_dFCHMM_rowSums_cpp,          1},
    {"_dFCHMM_standardize_rows_cpp", (DL_FUNC) &_dFCHMM_standardize_rows_cpp, 1},
    {NULL, NULL, 0}
};

void R_init_dFCHMM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
