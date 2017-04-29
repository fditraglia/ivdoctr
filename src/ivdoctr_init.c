#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ivdoctr_rinvwish(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ivdoctr_rinvwish", (DL_FUNC) &ivdoctr_rinvwish, 3},
    {NULL, NULL, 0}
};

void R_init_ivdoctr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
