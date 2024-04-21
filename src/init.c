
#include <R_ext/Visibility.h>
#include <Rinternals.h>
#include <Rversion.h>
#include "dgBIV.h"

// Use straightforward macro definitions for clarity and simplicity.
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(dgBIV, 4),
    {NULL, NULL, 0} // Terminating entry
};

/**
 * Initializes the genSurv package, registering routines and ensuring the
 * package namespace is correctly identified.
 *
 * @param dll Pointer to the DLL information.
 */
void attribute_visible R_init_genSurv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

    // Enforce symbol registration for R >= 2.16.0
    #if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
    #endif

    // Ensure the 'genSurv' namespace is correctly loaded and identified.
    SEXP genSurv_NS = R_FindNamespace(mkString("genSurv"));
    if (genSurv_NS == R_UnboundValue) {
        error("Missing 'genSurv' namespace: should never happen.");
    }
    if (!isEnvironment(genSurv_NS)) {
        error("'genSurv' namespace not determined correctly.");
    }
}
