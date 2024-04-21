validateGenCMMInputs <- function(n, model.cens, cens.par, beta, covar, rate) {
  if (n <= 0) stop("Argument 'n' must be greater than 0", call. = FALSE)
  if (!(model.cens %in% c("uniform", "exponential")))
    stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call. = FALSE)
  if (cens.par <= 0) stop("Argument 'cens.par' must be greater than 0", call. = FALSE)
  if (length(beta) != 3) stop("Argument 'beta' must be a vector with length 3", call. = FALSE)
  if (covar <= 0) stop("Argument 'covar' must be greater than 0", call. = FALSE)
  if (length(rate) != 6) stop("Argument 'rate' must be a vector with length 6", call. = FALSE)
}

generateEventTimes <- function(z1, beta, rate) {
  u <- runif(1, 0, 1)
  t12 <- (-log(1 - u) / (rate[1] * exp(beta[1] * z1)))^(1 / rate[2])
  u <- runif(1, 0, 1)
  t13 <- (-log(1 - u) / (rate[3] * exp(beta[2] * z1)))^(1 / rate[4])
  u <- runif(1, 0, 1)
  t23 <- (-log(1 - u) / (rate[5] * exp(beta[3] * z1)))^(1 / rate[6])
  
  return(list(t12 = t12, t13 = t13, t23 = t23))
}

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


genCMM <- function(n, model.cens, cens.par, beta, covar, rate) {
  validateGenCMMInputs(n, model.cens, cens.par, beta, covar, rate)
  
  data <- processEvents(n, cens.par, model.cens, beta, covar, rate)
  data <- data[-1, ]  # Remove the placeholder first row
  row.names(data) <- as.integer(1:nrow(data))
  class(data) <- c("data.frame", "CMM")
  
  return(data)
} # genCMM

