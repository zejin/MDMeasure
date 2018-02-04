#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h> // for NULL

/* declarations to register native routines in this package */ 

/* .C calls */
extern void est_complete(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"est_complete", (DL_FUNC) &est_complete, 7},
  {NULL, NULL, 0}
};

void R_init_energy(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}