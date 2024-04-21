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
#' validateGenTDCMInputs(n = 100, dist = "weibull", corr = 0.5, dist.par = c(1, 2, 3, 4),
#' model.cens = "uniform", cens.par = 0.5, beta = c(0.1, 0.2), lambda = 1)
validateGenTDCMInputs <- function(n, dist, corr, dist.par, model.cens, cens.par, beta, lambda) { # nolint
  if (n <= 0) stop("Argument 'n' must be greater than 0", call. = FALSE)
  if (!(dist %in% c("weibull", "exponential"))) 
    stop("Argument 'dist' must be one of 'weibull' or 'exponential'", call. = FALSE) # nolint
  if (dist == "weibull" && (corr <= 0 || corr > 1 || length(dist.par) != 4 || any(dist.par <= 0)))  # nolint
    stop("With 'dist=weibull', 'corr' must be in (0,1], and 'dist.par' must be a positive length-4 vector", call. = FALSE) # nolint
  if (dist == "exponential" && (corr < -1 || corr > 1 || length(dist.par) != 2 || any(dist.par <= 0))) # nolint
    stop("With dist='exponential', 'corr' must be in [-1,1], and 'dist.par' must be a positive length-2 vector", call. = FALSE) # nolint
  if (!(model.cens %in% c("uniform", "exponential"))) 
    stop("Argument 'model.cens' must be one of 'uniform' or 'exponential'", call. = FALSE) # nolint
  if (cens.par <= 0 || length(beta) != 2 || lambda <= 0) 
    stop("Arguments 'cens.par', 'beta', and 'lambda' must be positive, with 'beta' being a length-2 vector", call. = FALSE) # nolint
}

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

assembleDataFrame <- function(observations) {
  data <- data.frame(observations, row.names = NULL)
  names(data) <- c("id", "start", "stop", "event", "covariate", "tdcov")
  class(data) <- c("data.frame", "TDCM")
  return(data[-1, ])  # Exclude the placeholder first row
}


genTDCM <- function(n, dist, corr, dist.par, model.cens, cens.par, beta, lambda) {
  validateGenTDCMInputs(n, dist, corr, dist.par, model.cens, cens.par, beta, lambda)
  
  b <- dgBIV(n, dist, corr, dist.par)  # Assuming dgBIV is predefined or provided
  
  observations <- generateCensoredObservations(n, dist.par, model.cens, cens.par, beta, lambda, b)
  
  data <- assembleDataFrame(observations)
  
  row.names(data) <- as.integer(1:nrow(data))
  
  return(data)
} # genTDCM
