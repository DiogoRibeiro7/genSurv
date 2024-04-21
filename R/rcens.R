# Generate a vector of uniformly distributed censored values
# Args:
#   n: The number of observations to generate.
#   cens.par: The parameter to censor the uniform distribution at (max value).
# Returns:
#   A numeric vector of uniformly distributed values with specified censoring.
runifcens <- function(n, cens.par = 1) {
    # Ensure that cens.par is positive
    if (cens.par <= 0) stop("cens.par must be positive.")
    runif(n = n, min = 0, max = cens.par)
}

# Generate a vector of exponentially distributed censored values
# Args:
#   n: The number of observations to generate.
#   cens.par: The parameter to censor the exponential distribution at (rate parameter).
# Returns:
#   A numeric vector of exponentially distributed values with specified censoring.
rexpcens <- function(n, cens.par = 1) {
    # Ensure that cens.par is positive
    if (cens.par <= 0) stop("cens.par must be positive.")
    rexp(n = n, rate = 1 / cens.par)
}
