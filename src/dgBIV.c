
#include <R_ext/Random.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "defines.h"

// Define pointer types for clearer function signatures.
typedef const double* CdoubleCP;
typedef double* doubleCP;

// Typedef for the distribution function pointer for improved readability.
typedef void (*DistFunc)(CdoubleCP, CdoubleCP, doubleCP, doubleCP);

/**
 * Exponential distribution function.
 *
 * @param pcorr Pointer to the correlation coefficient.
 * @param pdistpar Pointer to the distribution parameters array.
 * @param t1 Pointer to the first return variable.
 * @param t2 Pointer to the second return variable.
 */
static void expt(CdoubleCP pcorr, CdoubleCP pdistpar, doubleCP t1, doubleCP t2) {
    double u1 = runif(0, 1);
    double v = runif(0, 1);
    double a = *pcorr * (2 * u1 - 1);
    double u2 = 2 * v / (1 - a + sqrt(pow(1 - a, 2) + 4 * a * v));
    *t1 = -pdistpar[0] * log(1 - u1);
    *t2 = -pdistpar[1] * log(1 - u2);
}

/**
 * Weibull distribution function.
 *
 * @param pcorr Pointer to the correlation coefficient.
 * @param pdistpar Pointer to the distribution parameters array.
 * @param t1 Pointer to the first return variable.
 * @param t2 Pointer to the second return variable.
 */
static void weibullt(CdoubleCP pcorr, CdoubleCP pdistpar, doubleCP t1, doubleCP t2) {
    double u[5];
    for (int i = 0; i < 5; i++) {
        u[i] = runif(0, 1);
    }
    double v = (u[4] > *pcorr) ? -log(u[3]) : (-log(u[1]) - log(u[2]));
    *t1 = pow(u[0], *pcorr / pdistpar[0]) * pow(v, 1 / pdistpar[0]) * pdistpar[1];
    *t2 = pow(1 - u[0], *pcorr / pdistpar[2]) * pow(v, 1 / pdistpar[2]) * pdistpar[3];
}

/**
 * Main function to generate a matrix of bivariate distributions.
 *
 * @param n Number of samples to generate.
 * @param dist Distribution type ("weibull" or "exponential").
 * @param corr Correlation coefficient.
 * @param distpar Distribution parameters.
 * @return SEXP Matrix of generated samples.
 */
SEXP dgBIV(SEXP n, SEXP dist, SEXP corr, SEXP distpar) {
    const char* pdist = CHAR(STRING_ELT(dist, 0));
    DistFunc tfunc = (strcmp(pdist, "weibull") == 0) ? weibullt : expt;

    SEXP mat;
    PROTECT(mat = allocMatrix(REALSXP, INTEGER(n)[0], 2));

    GetRNGstate();
    for (int i = 0; i < INTEGER(n)[0]; i++) {
        tfunc(REAL(corr), REAL(distpar), &REAL(mat)[i], &REAL(mat)[i + INTEGER(n)[0]]);
    }
    PutRNGstate();

    UNPROTECT(1);
    return mat;
}

