#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void LambdaC(void *, void *, void *, void *, void *, void *, void *, void *);
extern void LambdaK(void *, void *, void *, void *, void *, void *, void *, void *);
extern void Encapped(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sumpkC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"LambdaC",     (DL_FUNC) &LambdaC,      8},
    {"LambdaK",     (DL_FUNC) &LambdaK,      8},
    {"Encapped",     (DL_FUNC) &Encapped,    11},
    {"sumpkC",     (DL_FUNC) &sumpkC,      10},
    {NULL, NULL, 0}
};

void R_init_secrdesign(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
