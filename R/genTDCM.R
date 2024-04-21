#' Validate Input Parameters for TDCM Data Generation
#'
#' Ensures that the parameters provided for generating TDCM data are valid. It checks
#' the type, range, and length of each parameter according to the specifications for
#' the TDCM data generation process.
#'
#' @param n Integer; the number of observations to generate.
#' @param dist Character; specifies the distribution type, either "weibull" or "exponential".
#' @param corr Numeric; the correlation coefficient, constraints depend on `dist`.
#' @param dist.par Numeric vector; distribution parameters, their length and constraints depend on `dist`.
#' @param model.cens Character; the censoring model used, either "uniform" or "exponential".
#' @param cens.par Numeric; parameter for the censoring model, must be positive.
#' @param beta Numeric vector; of length 2, parameters influencing the model.
#' @param lambda Numeric; a positive parameter influencing the rate or intensity in the model.
#'
#' @return Does not return a value but stops function execution with an informative error
#' if any input is invalid.
#' @export
#' @examples
#' validateGenTDCMInputs(
#'   n = 100, dist = "weibull", corr = 0.5, dist.par = c(1, 2, 3, 4),
#'   model.cens = "uniform", cens.par = 0.5, beta = c(0.1, 0.2), lambda = 1
#' )
validateGenTDCMInputs <- function(n, dist, corr, dist.par, model.cens, cens.par, beta, lambda) { # nolint
  if (n <= 0) stop("Argument 'n' must be greater than 0", call. = FALSE)
  if (!(dist %in% c("weibull", "exponential"))) {
    stop("Argument 'dist' must be one of 'weibull' or 'exponential'", call. = FALSE)
  } # nolint
  if (dist == "weibull" && (corr <= 0 || corr > 1 || length(dist.par) != 4 || any(dist.par <= 0))) { # nolint
    stop("With 'dist=weibull', 'corr' must be in (0,1], and 'dist.par' must be a positive length-4 vector", call. = FALSE)
  } # nolint
  if (dist == "exponential" && (corr < -1 || corr > 1 || length(dist.par) != 2 || any(dist.par <= 0))) { # nolint
    stop("With dist='exponential', 'corr' must be in [-1,1], and 'dist.par' must be a positive length-2 vector", call. = FALSE)
  } # nolint
  if (!(model.cens %in% c("uniform", "exponential"))) {
    stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call. = FALSE)
  } # nolint
  if (cens.par <= 0 || length(beta) != 2 || lambda <= 0) {
    stop("Arguments 'cens.par', 'beta', and 'lambda' must be positive, with 'beta' being a length-2 vector", call. = FALSE)
  } # nolint
}

#' Generate Censored Observations Based on Specified Parameters
#'
#' Simulates censored observations for a dataset based on the specified distribution,
#' censoring model, and other parameters. It is designed to support both uniform and
#' exponential censoring models.
#'
#' @param n Integer; the number of observations to generate.
#' @param dist.par Numeric vector; parameters of the underlying distribution for the
#'        time to event data. The interpretation and necessary length of this vector
#'        depend on the distribution being simulated.
#' @param model.cens Character; specifies the censoring model to use, either "uniform"
#'        or "exponential".
#' @param cens.par Numeric; parameter for the censoring model, must be positive.
#' @param beta Numeric vector of length 2; parameters influencing the rate or intensity
#'        in the model.
#' @param lambda Numeric; a positive parameter influencing the rate or intensity in the model.
#' @param b Matrix; a matrix of covariate values where each row corresponds to an observation
#'        and the first two columns are used in the simulation process.
#'
#' @return A matrix where each row represents a censored observation with columns for:
#'         observation index, start time (always 0 in this context), end time (censored or event time),
#'         event status (1 if the event occurred, 0 if censored), and two covariate values (z1 and z2).
#'         The column names are "id", "start", "stop", "status", "covariate1", "covariate2".
#' @examples
#' n <- 100
#' dist.par <- c(1, 2) # Example parameters
#' model.cens <- "uniform"
#' cens.par <- 0.5
#' beta <- c(0.1, 0.2)
#' lambda <- 1
#' b <- matrix(rnorm(n * 2), ncol = 2)
#' observations <- generateCensoredObservations(n, dist.par, model.cens, cens.par, beta, lambda, b)
#' @export
generateCensoredObservations <- function(n, dist.par, model.cens, cens.par, beta, lambda, b) {
  rfunc <- if (model.cens == "uniform") runifcens else rexpcens
  observations <- vector("list", n)

  for (k in 1:n) {
    z1 <- b[k, 2]
    c <- rfunc(1, cens.par)
    u <- runif(1)
    if (u < 1 - exp(-lambda * b[k, 1] * exp(beta[1] * z1))) {
      t <- -log(1 - u) / (lambda * exp(beta[1] * z1))
      z2 <- 0
    } else {
      t <- -(log(1 - u) + lambda * b[k, 1] * exp(beta[1] * z1) * (1 - exp(beta[2]))) / (lambda * exp(beta[1] * z1 + beta[2]))
      z2 <- 1
    }
    time <- min(t, c)
    status <- ifelse(t > c, 0, 1)
    observations[[k]] <- list(c(k, 0, time, status, z1, z2))
  }

  return(do.call(rbind, observations))
}

#' Assemble Observations into a Data Frame with Specified Structure
#'
#' Takes a list of observations and converts it into a structured data frame,
#' applying the appropriate class and column names. It is specifically tailored
#' for assembling Time-Dependent Covariate Model (TDCM) data frames from observation lists.
#'
#' @param observations List; a list of observations where each element is expected
#'        to be a numeric vector with the elements representing id, start time,
#'        stop time, event occurrence, and covariates. The list format is typically
#'        generated from simulation functions.
#'
#' @return A data frame of the assembled observations with columns named "id",
#'         "start", "stop", "event", "covariate", and "tdcov". The data frame is
#'         also assigned the class "TDCM" in addition to "data.frame". The first row,
#'         often used as a placeholder in simulation functions, is excluded from the
#'         returned data frame.
#' @examples
#' observations <- list(
#'   c(1, 0, 2.5, 1, 0.5, 0),
#'   c(2, 0, 3.0, 0, 0.7, 1)
#' )
#' df <- assembleDataFrame(observations)
#' print(df)
#' @export
assembleDataFrame <- function(observations) {
  data <- data.frame(observations, row.names = NULL)
  names(data) <- c("id", "start", "stop", "event", "covariate", "tdcov")
  class(data) <- c("data.frame", "TDCM")
  return(data[-1, ]) # Exclude the placeholder first row
}

#' Generate Time-Dependent Covariate Model (TDCM) Data
#'
#' Simulates data for a Time-Dependent Covariate Model (TDCM) based on the specified
#' distribution for event times, correlation parameters, censoring model, and other
#' parameters affecting the simulation process.
#'
#' @param n Integer; the number of observations to generate.
#' @param dist Character; the distribution of event times, either "weibull" or "exponential".
#' @param corr Numeric; the correlation parameter, relevant for the distribution of event times.
#' @param dist.par Numeric vector; parameters of the specified distribution for event times.
#' @param model.cens Character; the censoring model, either "uniform" or "exponential".
#' @param cens.par Numeric; the parameter for the censoring model, must be positive.
#' @param beta Numeric vector; parameters influencing the intensity function in the model.
#' @param lambda Numeric; the rate parameter for the intensity function, must be positive.
#'
#' @return A `data.frame` of class "TDCM" containing the generated observations.
#'         Each row corresponds to an observation with columns for the observation id,
#'         start and stop times, event occurrence indicator, and covariates. 
#'
#' @examples
#' n <- 100
#' dist <- "weibull"
#' corr <- 0.5
#' dist.par <- c(2, 2, 2, 2) # Example for Weibull distribution
#' model.cens <- "uniform"
#' cens.par <- 0.5
#' beta <- c(0.1, 0.2)
#' lambda <- 1
#' data <- genTDCM(n, dist, corr, dist.par, model.cens, cens.par, beta, lambda)
#'
#' @export
genTDCM <- function(n, dist, corr, dist.par, model.cens, cens.par, beta, lambda) {
  validateGenTDCMInputs(n, dist, corr, dist.par, model.cens, cens.par, beta, lambda)

  b <- dgBIV(n, dist, corr, dist.par) # Assuming dgBIV is predefined or provided

  observations <- generateCensoredObservations(n, dist.par, model.cens, cens.par, beta, lambda, b)

  data <- assembleDataFrame(observations)

  row.names(data) <- as.integer(1:nrow(data))

  return(data)
} # genTDCM
