#' Validate Input Parameters for CPHM Data Generation
#'
#' Ensures that the provided parameters meet the required criteria for generating
#' Cox Proportional Hazards Model (CPHM) data with censoring.
#'
#' @param n Integer; the number of data points to generate.
#' @param model.cens Character; specifies the censoring model, either "uniform" or "exponential".
#' @param cens.par Numeric; the parameter for the censoring model, must be positive.
#' @param covar Numeric; the covariate value, must be positive.
#'
#' @return Does not return a value but stops execution with an error message if any input is invalid.
#' @export
#' @examples
#' validateGenCPHMInputs(n = 100, model.cens = "uniform", cens.par = 1, covar = 0.5)
validateGenCPHMInputs <- function(n, model.cens, cens.par, covar) {
    if (n <= 0) stop("Argument 'n' must be greater than 0", call. = FALSE)
    if (!(model.cens %in% c("uniform", "exponential"))) {
        stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call. = FALSE)
    }
    if (cens.par <= 0) stop("Argument 'cens.par' must be greater than 0", call. = FALSE)
    if (covar <= 0) stop("Argument 'covar' must be greater than 0", call. = FALSE)
}

#' Generate CPHM Data
#'
#' Simulates data for the Cox Proportional Hazards Model (CPHM) based on specified parameters
#' and a censoring model.
#'
#' @param n Integer; the number of observations to generate.
#' @param rfunc Function; a function to generate censored times, determined by `model.cens`.
#' @param cens.par Numeric; parameter for the censoring model.
#' @param beta Numeric; the coefficient for the covariate in the model.
#' @param covar Numeric; the range of the covariate values.
#'
#' @return A matrix with `n` rows and 3 columns: time to event or censoring, event indicator (1 for event, 0 for censoring), and covariate value.
#' @export
#' @examples
#' n <- 100
#' beta <- 0.5
#' covar <- 2
#' data <- generateCPHMData(n, runifcens, 1, beta, covar)
#' data <- data.frame(data)
#' names(data) <- c("time", "status", "covariate")
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

#' Generate Data Using the Cox Proportional Hazards Model (CPHM)
#'
#' This function generates simulated data based on the Cox Proportional Hazards Model (CPHM),
#' incorporating a specified censoring model and covariate effects.
#'
#' @param n Integer; the number of observations to generate.
#' @param model.cens Character; the censoring model, either "uniform" or "exponential".
#' @param cens.par Numeric; parameter for the censoring model.
#' @param beta Numeric; the coefficient for the covariate in the model.
#' @param covar Numeric; the range of the covariate values.
#'
#' @return A `data.frame` marked with class "CPHM" containing the generated observations.
#'         Each row corresponds to an observation with columns for time to event or censoring,
#'         event status, and the covariate value.
#' @export
#' @examples
#' n <- 100
#' model.cens <- "uniform"
#' cens.par <- 1
#' beta <- 0.5
#' covar <- 2
#' data <- genCPHM(n, model.cens, cens.par, beta, covar)
#' print(head(data))
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
