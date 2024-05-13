#' Validates the inputs for the genCMM function.
#'
#' This function validates the inputs for the genCMM function, which generates
#' survival data based on a continuous-time multi-state Markov model.
#'
#' @param n The number of individuals in the generated data.
#' @param model.cens The censoring model to be used. Must be one of "uniform" or "exponential".
#' @param cens.par The parameter for the censoring model. Must be greater than 0.
#' @param beta The vector of regression coefficients for the covariates. Must have length 3.
#' @param covar The covariate value. Must be greater than 0.
#' @param rate The vector of transition rates between states. Must have length 6.
#'
#' @return None
#'
#' @examples
#' validateGenCMMInputs(n = 100, model.cens = "uniform", cens.par = 1, beta = c(0.5, 0.2, -0.1),
#'                      covar = 2, rate = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
#'
#' @export
validateGenCMMInputs <- function(n, model.cens, cens.par, beta, covar, rate) {
  if (n <= 0) stop("Argument 'n' must be greater than 0", call. = FALSE)
  if (!(model.cens %in% c("uniform", "exponential")))
    stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call. = FALSE)
  if (cens.par <= 0) stop("Argument 'cens.par' must be greater than 0", call. = FALSE)
  if (length(beta) != 3) stop("Argument 'beta' must be a vector with length 3", call. = FALSE)
  if (covar <= 0) stop("Argument 'covar' must be greater than 0", call. = FALSE)
  if (length(rate) != 6) stop("Argument 'rate' must be a vector with length 6", call. = FALSE)
}

# Function to generate event times based on a set of parameters
# 
# Args:
#   z1: Numeric value representing the input variable z1
#   beta: Numeric vector representing the beta parameters
#   rate: Numeric vector representing the rate parameters
# 
# Returns:
#   A list containing the generated event times t12, t13, and t23
generateEventTimes <- function(z1, beta, rate) {
  u <- runif(1, 0, 1)
  t12 <- (-log(1 - u) / (rate[1] * exp(beta[1] * z1)))^(1 / rate[2])
  u <- runif(1, 0, 1)
  t13 <- (-log(1 - u) / (rate[3] * exp(beta[2] * z1)))^(1 / rate[4])
  u <- runif(1, 0, 1)
  t23 <- (-log(1 - u) / (rate[5] * exp(beta[3] * z1)))^(1 / rate[6])
  
  return(list(t12 = t12, t13 = t13, t23 = t23))
}

#' Generate Survival Data with Competing Risks
#'
#' This function generates survival data with competing risks based on the specified parameters.
#'
#' @param n The number of individuals to generate data for.
#' @param cens.par The parameter for the censoring distribution.
#' @param model.cens The type of censoring distribution ("uniform" or "exponential").
#' @param beta The coefficient for the covariate in the event time generation model.
#' @param covar The range of the covariate values.
#' @param rate The rate parameter for the event time generation model.
#'
#' @return A data frame containing the generated survival data with columns for individual ID, start time, stop time,
#' event type, covariate value, and transition type.
#'
#' @examples
#' data <- processEvents(100, 1, "uniform", 0.5, 10, 0.1)
#' head(data)
#'
#' @export
processEvents <- function(n, cens.par, model.cens, beta, covar, rate) {
  rfunc <- if (model.cens == "uniform") runifcens else rexpcens
  mat <- list()  # Use a list to collect rows; more efficient than rbind in a loop

  for (k in 1:n) {
    z1 <- runif(1, 0, covar)
    c <- rfunc(1, cens.par)
    events <- generateEventTimes(z1, beta, rate)
    t12 <- events$t12
    t13 <- events$t13
    t23 <- events$t23

    # Determine the earliest event time and type
    minEventTime <- min(t12, t13, c)
    eventType <- ifelse(minEventTime == t12, 1,
           ifelse(minEventTime == t13, 2, NA))
    eventObserved <- minEventTime < c

    # Record the event or censoring
    if (!is.na(eventType) && eventObserved) {
      # Event is observed before censoring
      mat[[k]] <- c(k, 0, minEventTime, eventType, z1, eventType)
    } else {
      # Censoring occurs before any event or no event type is determined
      eventType <- 0  # Indicating censoring
      mat[[k]] <- c(k, 0, c, eventType, z1, NA)  # NA for transition type in case of censoring
    }
  }

  # Convert the list to a matrix and then to a dataframe
  mat <- do.call(rbind, mat)
  data <- data.frame(mat, stringsAsFactors = FALSE)
  names(data) <- c("id", "start", "stop", "event", "covariate", "trans")
  row.names(data) <- NULL
  return(data)
}


#' Generate a Competing Risks Mixture Model (CMM) dataset
#'
#' This function generates a Competing Risks Mixture Model (CMM) dataset based on the provided parameters.
#'
#' @param n The number of observations to generate.
#' @param model.cens The censoring model to use. Can be "exponential" or "weibull".
#' @param cens.par The parameters for the censoring model. For "exponential", it is the rate parameter. For "weibull", it is a vector of shape and scale parameters.
#' @param beta The regression coefficients for the covariates.
#' @param covar The covariate matrix.
#' @param rate The baseline hazard rate.
#'
#' @return A data frame representing the generated CMM dataset.
#'
#' @examples
#' genCMM(n = 100, model.cens = "exponential", cens.par = 0.1, beta = c(1, 2), covar = matrix(c(1, 2, 3, 4), ncol = 2), rate = 0.05)
#'
#' @export
genCMM <- function(n, model.cens, cens.par, beta, covar, rate) {
  validateGenCMMInputs(n, model.cens, cens.par, beta, covar, rate)
  
  data <- processEvents(n, cens.par, model.cens, beta, covar, rate)
  data <- data[-1, ]  # Remove the placeholder first row
  row.names(data) <- as.integer(1:nrow(data))
  class(data) <- c("data.frame", "CMM")
  
  return(data)
} # genCMM

