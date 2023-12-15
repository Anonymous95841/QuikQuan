#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/RS.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/
#define FDEF(name){#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}

/* .Fortran calls */
void F77_NAME(fast_kqr)(double *delta, double *Kmat, int *nobs,
    double *y, int *nlam, double *ulam, double *tau, double *eps, int *maxit,
    double *gamma, int *is_exact, int *anlam, int *npass, int *jerr,
    double *alpmat);

static R_NativePrimitiveArgType fast_kqr_t[] = {REALSXP, REALSXP, INTSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP};

void F77_NAME(fast_nckqr)(double *delta, double *Kmat, int *nobs,
    double *y, int *nlam1, double *ulam1, int *nlam2, double *ulam2,
    int *ntau, double *utau, double *eps, int *maxit, double *gamma,
    int *is_exact, int *anlam1, int *anlam2, int *npass, int *jerr,
    double *alpmat);

static R_NativePrimitiveArgType fast_nckqr_t[] = {REALSXP, REALSXP, INTSXP,
  REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP,
  INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP};

void F77_NAME(rbfdot)(double *X1, double *X2, int *nobs1, int *nobs2,
    int *p1, int *p2, double *sigma, double *Kmat);

static R_NativePrimitiveArgType rbfdot_t[] = {REALSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static R_FortranMethodDef FortranEntries[] = {
    FDEF(fast_kqr) ,
    FDEF(fast_nckqr) ,
    FDEF(rbfdot),
    {NULL, NULL, 0}
};

void R_init_QuikQuan(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
