#' Validate the inputs for the dgBIV function
#'
#' This function validates the inputs for the dgBIV function, which generates
#' random samples from a bivariate distribution.
#'
#' @param n The number of samples to generate. Must be a positive integer.
#' @param dist The distribution to use. Must be one of "weibull" or "exponential".
#' @param corr The correlation coefficient between the two variables. Must be a numeric value between -1 and 1.
#' @param dist.par A numeric vector of parameters for the distribution. Must contain positive values.
#'
#' @return This function does not return anything. It throws an error if any of the inputs are invalid.
#'
#' @examples
#' validateDgBIVInputs(100, "weibull", 0.5, c(1, 2))
#' validateDgBIVInputs(200, "exponential", -0.8, c(3, 4, 5))
#' validateDgBIVInputs(500, "normal", 0.2, c(2))
#' validateDgBIVInputs(1000, "weibull", 0.9, c(-1, 2, 3))
#' validateDgBIVInputs(100, "exponential", 0.5, c(1, 2, 3, 0))
#'
#' @seealso \code{\link{dgBIV}}
#' @keywords validation
validateDgBIVInputs <- function(n, dist, corr, dist.par) {
    if (!is.numeric(n) || length(n) != 1 || n <= 0) {
        stop("Argument 'n' must be a positive integer", call. = FALSE)
    }
    if (!dist %in% c("weibull", "exponential")) {
        stop("Argument 'dist' must be one of 'weibull' or 'exponential'", call. = FALSE)
    }
    if (!is.numeric(corr) || length(corr) != 1 || corr <= -1 || corr >= 1) {
        stop("Argument 'corr' must be a numeric value between -1 and 1", call. = FALSE)
    }
    if (!is.numeric(dist.par) || length(dist.par) < 1) {
        stop("Argument 'dist.par' must be a numeric vector with positive values", call. = FALSE)
    }
    if (any(dist.par <= 0)) {
        stop("All elements in 'dist.par' must be greater than 0", call. = FALSE)
    }
}

#' Calculate the density of a bivariate distribution
#'
#' This function calculates the density of a bivariate distribution given the number of samples, the distribution type, the correlation coefficient, and the distribution parameters.
#'
#' @param n The number of samples to generate.
#' @param dist The type of distribution. This should be a string specifying the distribution type (e.g., "normal", "uniform").
#' @param corr The correlation coefficient between the two variables. This should be a double value between -1 and 1.
#' @param dist.par A vector of distribution parameters. The length and meaning of the parameters depend on the distribution type.
#'
#' @return The density of the bivariate distribution.
#'
#' @examples
#' dgBIV(100, "normal", 0.5, c(0, 1))
#'
#' @export
dgBIV <- function(n, dist, corr, dist.par) {
    validateDgBIVInputs(n, dist, corr, dist.par)

    # n (integer), dist (string), corr (double), dist.par (double vector)
    .Call(Rf_dgBIV, as.integer(n), dist, as.double(corr), as.double(dist.par), PACKAGE="genSurv")
}
