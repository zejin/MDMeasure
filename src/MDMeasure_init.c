#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h> // for NULL

////////////////////////////////////////////////////////////////////////////////////////////////////
///// declarations to register native routines in this package
////////////////////////////////////////////////////////////////////////////////////////////////////

/* .C calls */
extern void dCov(double *X, double *Y, double *XX, double *YY, double *Q, int *NOBS, int *NDIMX, int *NDIMY);
extern void dCov_perm(double *XX, double *YY, double *Q, int *NOBS, int *IPERM);

extern void dCov_asymmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
extern void dCov_asymmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);

extern void dCov_symmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
extern void dCov_symmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);

extern void MDM_complete(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
extern void MDM_complete_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);

extern void MDM_complete_simple(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
extern void MDM_complete_simple_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);

extern void MDM_asymmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
extern void MDM_asymmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);

extern void MDM_asymmetric_simple(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
extern void MDM_asymmetric_simple_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);

extern void MDM_symmetric(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
extern void MDM_symmetric_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);

extern void MDM_symmetric_simple(double *X, double *D, double *Q, int *NOBS, int *NDIM, int *NCOMP, int *ICOMP);
extern void MDM_symmetric_simple_perm(double *D, double *Q, int *NOBS, int *NCOMP, int *IPERM);

static const R_CMethodDef CEntries[] = {
  {"dCov",                       (DL_FUNC) &dCov,                       8},
  {"dCov_perm",                  (DL_FUNC) &dCov_perm,                  5},
  {"dCov_asymmetric",            (DL_FUNC) &dCov_asymmetric,            7},
  {"dCov_asymmetric_perm",       (DL_FUNC) &dCov_asymmetric_perm,       5},
  {"dCov_symmetric",             (DL_FUNC) &dCov_symmetric,             7},
  {"dCov_symmetric_perm",        (DL_FUNC) &dCov_symmetric_perm,        5},
  {"MDM_complete",               (DL_FUNC) &MDM_complete,               7},
  {"MDM_complete_perm",          (DL_FUNC) &MDM_complete_perm,          5},
  {"MDM_complete_simple",        (DL_FUNC) &MDM_complete_simple,        7},
  {"MDM_complete_simple_perm",   (DL_FUNC) &MDM_complete_simple_perm,   5},
  {"MDM_asymmetric",             (DL_FUNC) &MDM_asymmetric,             7},
  {"MDM_asymmetric_perm",        (DL_FUNC) &MDM_asymmetric_perm,        5},
  {"MDM_asymmetric_simple",      (DL_FUNC) &MDM_asymmetric_simple,      7},
  {"MDM_asymmetric_simple_perm", (DL_FUNC) &MDM_asymmetric_simple_perm, 5},
  {"MDM_symmetric",              (DL_FUNC) &MDM_symmetric,              7},
  {"MDM_symmetric_perm",         (DL_FUNC) &MDM_symmetric_perm,         5},
  {"MDM_symmetric_simple",       (DL_FUNC) &MDM_symmetric_simple,       7},
  {"MDM_symmetric_simple_perm",  (DL_FUNC) &MDM_symmetric_simple_perm,  5},
  {NULL, NULL, 0}
};

void R_init_energy(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}