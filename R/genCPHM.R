validateGenCPHMInputs <- function(n, model.cens, cens.par, covar) {
    if (n <= 0) stop("Argument 'n' must be greater than 0", call. = FALSE)
    if (!(model.cens %in% c("uniform", "exponential"))) {
        stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call. = FALSE)
    }
    if (cens.par <= 0) stop("Argument 'cens.par' must be greater than 0", call. = FALSE)
    if (covar <= 0) stop("Argument 'covar' must be greater than 0", call. = FALSE)
}

generateCPHMData <- function(n, rfunc, cens.par, beta, covar) {
    data <- matrix(data = 0, ncol = 3, nrow = n)

    for (k in 1:n) {
        z <- runif(1, 0, covar)
        c <- rfunc(1, cens.par)
        x <- rexp(1, exp(beta * z))

        time <- min(x, c)
        status <- ifelse(x > c, 0, 1)

        data[k, ] <- c(time, status, z)
    }

    return(data)
}

genCPHM <- function(n, model.cens, cens.par, beta, covar) {
    validateGenCPHMInputs(n, model.cens, cens.par, covar)

    rfunc <- if (model.cens == "uniform") runifcens else rexpcens
    data <- generateCPHMData(n, rfunc, cens.par, beta, covar)

    data <- data.frame(data)
    names(data) <- c("time", "status", "covariate")
    row.names(data) <- as.integer(1:nrow(data))
    class(data) <- c("data.frame", "CPHM")

    return(data)
} # genCPHM
