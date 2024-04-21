#' Validate input parameters for THMM data generation
#'
#' Checks if the provided arguments meet the expected criteria for generating
#' time-homogeneous Markov model (THMM) data.
#'
#' @param n Integer; the number of data points to generate.
#' @param model.cens Character; specifies the censoring model, either "uniform" or "exponential".
#' @param cens.par Numeric; the parameter for the censoring model, must be positive.
#' @param beta Numeric vector of length 3; parameters for the beta distribution.
#' @param covar Numeric; the covariate value, must be positive.
#' @param rate Numeric vector of length 3; transition rates for the Markov model.
#'
#' @return Invisibly returns NULL but will stop execution with an error message if inputs are invalid.
#' @export
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

#' Generate censored time based on specified model
#'
#' This function generates a single censored time value based on the specified
#' censoring model and its parameter.
#'
#' @param model.cens Character; the censoring model used, either "uniform" or "exponential".
#' @param cens.par Numeric; the parameter for the censoring model.
#'
#' @return Numeric; a single censored time value.
#' @export
getCensoredTime <- function(model.cens, cens.par) {
    if (model.cens == "uniform") {
        return(runifcens(1, cens.par))
    } else { # Exponential
        return(rexpcens(1, cens.par))
    }
}

#' Calculate transition times for THMM
#'
#' Computes the transition times between states in the time-homogeneous Markov
#' model based on given rates and covariate values.
#'
#' @param z1 Numeric; covariate value.
#' @param cens.par Numeric; censoring parameter.
#' @param beta Numeric vector; parameters influencing transition rates.
#' @param rate Numeric vector; base rates for transitions.
#' @param rfunc Function; function to generate censored time.
#'
#' @return List containing censored time 'c' and transition times 't12', 't13', 't23'.
#' @export
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

#' Assemble THMM data from individual records
#'
#' Constructs a data frame from a list of individual transition records, assigning
#' appropriate column names and class.
#'
#' @param records List; a list where each element is a vector representing a transition record.
#'
#' @return Data frame with THMM data, marked with class 'THMM'.
#' @export
assembleData <- function(records) {
    data <- data.frame(do.call(rbind, records), row.names = NULL)
    names(data) <- c("PTNUM", "time", "state", "covariate")
    data <- data[-1, ] # Remove initial placeholder row
    row.names(data) <- as.integer(1:nrow(data))
    class(data) <- c("data.frame", "THMM")
    return(data)
}

#' Generate THMM data
#'
#' Simulates data for a time-homogeneous Markov model based on specified parameters
#' and censoring model.
#'
#' @param n Integer; number of observations to generate.
#' @param model.cens Character; censoring model, either "uniform" or "exponential".
#' @param cens.par Numeric; parameter for the censoring model.
#' @param beta Numeric vector; influences transition rates.
#' @param covar Numeric; covariate value.
#' @param rate Numeric vector; base transition rates.
#'
#' @return Data frame of THMM data with columns for patient number, time, state, and covariate.
#' @export
#' @examples
#' genTHMM(n = 100, model.cens = "uniform", cens.par = 0.5, beta = c(0.1, 0.2, 0.3), covar = 1, rate = c(0.5, 0.6, 0.7))
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
} 
