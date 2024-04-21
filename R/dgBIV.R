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

dgBIV <- function(n, dist, corr, dist.par) {
    validateDgBIVInputs(n, dist, corr, dist.par)

    # n (integer), dist (string), corr (double), dist.par (double vector)
    .Call(Rf_dgBIV, as.integer(n), dist, as.double(corr), as.double(dist.par), PACKAGE="genSurv")
}

