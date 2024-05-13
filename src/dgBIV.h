
/**
 * Generate random samples from a bivariate distribution.
 *
 * This function generates random samples from a bivariate distribution specified by the given parameters.
 *
 * @param n       The number of samples to generate.
 * @param dist    The distribution type of the bivariate distribution.
 * @param corr    The correlation coefficient of the bivariate distribution.
 * @param distpar The distribution parameters of the bivariate distribution.
 *
 * @return        A SEXP object containing the generated random samples.
 */
SEXP dgBIV(SEXP n, SEXP dist, SEXP corr, SEXP distpar);
