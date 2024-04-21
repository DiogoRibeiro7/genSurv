validateInputs <- function(n, model.cens, cens.par, beta, covar, rate) { # nolint
    if (n <= 0) stop("Argument 'n' must be greater than 0", call. = FALSE)
    if (!(model.cens %in% c("uniform", "exponential"))) {
        stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call. = FALSE)
    } # nolint
    if (cens.par <= 0) stop("Argument 'cens.par' must be greater than 0", call. = FALSE) # nolint
    if (length(beta) != 3) stop("Argument 'beta' length must be a vector with length 3", call. = FALSE) # nolint
    if (covar <= 0) stop("Argument 'covar' must be greater than 0", call. = FALSE)
    if (length(rate) != 3) stop("Argument 'rate' must be a vector with length 3", call. = FALSE) # nolint
}


getCensoredTime <- function(model.cens, cens.par) {
    if (model.cens == "uniform") {
        return(runifcens(1, cens.par))
    } else { # Exponential
        return(rexpcens(1, cens.par))
    }
}

calculateTransitions <- function(z1, cens.par, beta, rate, rfunc) {
    c <- rfunc(1, cens.par)
    rate12 <- rate[1] * exp(beta[1] * z1)
    rate13 <- rate[2] * exp(beta[2] * z1)
    rate23 <- rate[3] * exp(beta[3] * z1)

    t12 <- rexp(1, rate12)
    t13 <- rexp(1, rate13)
    t23 <- rexp(1, rate23)

    return(list(c = c, t12 = t12, t13 = t13, t23 = t23))
}

assembleData <- function(records) {
    data <- data.frame(do.call(rbind, records), row.names = NULL)
    names(data) <- c("PTNUM", "time", "state", "covariate")
    data <- data[-1, ] # Remove initial placeholder row
    row.names(data) <- as.integer(1:nrow(data))
    class(data) <- c("data.frame", "THMM")
    return(data)
}

genTHMM <- function(n, model.cens, cens.par, beta, covar, rate) {
    validateInputs(n, model.cens, cens.par, beta, covar, rate)
    rfunc <- if (model.cens == "uniform") runifcens else rexpcens
    records <- vector("list", n)

    for (k in 1:n) {
        z1 <- runif(1, 0, covar)
        trans <- calculateTransitions(z1, cens.par, beta, rate, rfunc)
        c <- trans$c
        t12 <- trans$t12
        t13 <- trans$t13
        t23 <- trans$t23

        # Determine state based on transition times and censoring
        for (k in 1:n) {
            z1 <- runif(1, 0, covar)
            trans <- calculateTransitions(z1, cens.par, beta, rate, rfunc)
            c <- trans$c
            t12 <- trans$t12
            t13 <- trans$t13
            t23 <- trans$t23

            if (c < min(t12, t13)) {
                # Censoring occurs before any transition
                records[[k]] <- list(c(k, 0, 1, z1), c(k, c, 1, z1))
            } else {
                if (t13 < t12) {
                    # Transition to state 3 before state 2
                    records[[k]] <- list(c(k, 0, 1, z1), c(k, t13, 3, z1))
                } else {
                    # Transition to state 2 occurs first
                    records[[k]] <- list(c(k, 0, 1, z1), c(k, t12, 2, z1))
                    if (c < t12 + t23) {
                        # Censoring after transition to state 2 but before transition to state 3
                        records[[k]] <- c(records[[k]], list(c(k, c, 2, z1)))
                    } else {
                        # Transition to state 3 occurs after state 2 and is not censored
                        records[[k]] <- c(records[[k]], list(c(k, t12 + t23, 3, z1)))
                    }
                }
            }
        }
    }

    return(assembleData(records))
} # genTHMM
