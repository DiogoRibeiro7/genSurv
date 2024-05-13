#' Check if an object is of class TDCM
#'
#' This function checks if an object is a data frame and inherits the class "TDCM".
#'
#' @param x The object to be checked.
#' @return A logical value indicating whether the object is of class TDCM.
#' @examples
#' is.TDCM(data.frame()) # FALSE
#' is.TDCM(TDCM()) # TRUE
#' @export
is.TDCM <- function(x) {
    is.data.frame(x) && inherits(x, "TDCM")
}


#' Convert an object to TDCM class
#'
#' This function is used to convert an object to the TDCM (Time-Dependent Covariate Model) class.
#' It dispatches to the appropriate method based on the class of the input object.
#'
#' @param x The object to be converted.
#'
#' @return An object of class TDCM.
#'
#' @export
as.TDCM <- function(x) {
    UseMethod("as.TDCM")
}


#' Coerce an object into class 'TDCM'
#'
#' This function attempts to coerce an object into the class 'TDCM'. If the coercion is not possible, an error is thrown.
#'
#' @param x The object to be coerced.
#' @return The coerced object.
#' @export
as.TDCM.default <- function(x) {
    stop("cannot coerce class '", deparse(substitute(x)), "' into class 'TDCM'", domain = NA)
}


#' Convert to TDCM object
#'
#' This function converts an object to the TDCM class.
#'
#' @param x The object to be converted.
#' @return The converted object of class 'TDCM'.
#' @export
#'
#' @examples
#' as.TDCM.TDCM(1)
#'
#' @importFrom methods is
#' @importFrom methods stop
#' @importFrom methods is.TDCM
as.TDCM.TDCM <- function(x) {
    if (!is.TDCM(x)) stop("'x' must be of class 'TDCM'")
    x
}

#' Convert CMM object to TDCM object
#'
#' This function converts a CMM (Continuous-time Markov Model) object to a TDCM (Time-Dependent Covariate Model) object.
#'
#' @param x A CMM object.
#'
#' @return A TDCM object.
#'
#' @details The function takes a CMM object as input and returns a TDCM object. The TDCM object is created by transforming the rows of the input CMM object based on certain conditions. The transformed rows are then combined into a data frame and assigned the class "TDCM".
#'
#' @examples
#' # Create a CMM object
#' cmm <- create.CMM(...)
#'
#' # Convert CMM object to TDCM object
#' tdcmm <- as.TDCM(cmm)
#'
#' @export
as.TDCM.CMM <- function(x) {
    if (!is.CMM(x)) stop("'x' must be of class 'CMM'")

    # Preparing the initial data frame more directly
    data <- x[c("id", "start", "stop", "event", "covariate", "trans")]

    # Initialize an empty list for dynamic row addition
    transformedRows <- list()

    # Iterate over rows of data to transform
    for (i in 1:nrow(data)) {
        currentID <- data$id[i]
        covariate <- data$covariate[i]
        trans <- data$trans[i]
        event <- data$event[i]
        startTime <- ifelse(i == 1 || data$id[i] != data$id[i - 1], 0, data$stop[i - 1])
        stopTime <- data$stop[i]

        # Create transformed row(s) based on conditions
        if (event == 0 && (i == nrow(data) || data$id[i] != data$id[i + 1])) {
            transformedRows[[length(transformedRows) + 1]] <- c(currentID, startTime, stopTime, event, covariate, trans)
        } else if (event == 1) {
            transformedRows[[length(transformedRows) + 1]] <- c(currentID, startTime, stopTime, event, covariate, trans)
            if (i < nrow(data) && data$id[i] == data$id[i + 1]) {
                nextStopTime <- data$stop[i + 1]
                nextEvent <- data$event[i + 1]
                transformedRows[[length(transformedRows) + 1]] <- c(currentID, stopTime, nextStopTime, nextEvent, covariate, trans)
                i <- i + 1 # Skip next row as it's already processed
            }
        }
    }

    # Convert list to data frame
    data2 <- do.call(rbind, transformedRows)
    colnames(data2) <- c("id", "start", "stop", "event", "covariate", "tdcov")
    data2 <- data.frame(data2, stringsAsFactors = FALSE)
    class(data2) <- "TDCM"

    return(data2)
}


#' Convert a THMM object to a TDCM object
#'
#' This function takes a THMM (Time-Homogeneous Markov Model) object and converts it to a TDCM (Time-Dependent Covariate Model) object.
#'
#' @param x A THMM object.
#'
#' @return A TDCM object.
#'
#' @details The function transforms the data in the THMM object to match the format required for a TDCM object. It extracts the necessary columns from the THMM object, adjusts the rows based on state transitions, and converts the resulting list to a data frame with appropriate column names. The resulting data frame is then assigned the class "TDCM" and returned.
#'
#' @examples
#' # Create a THMM object
#' thmm <- create.THMM(...)
#'
#' # Convert the THMM object to a TDCM object
#' tdcmm <- as.TDCM.THMM(thmm)
#'
#' @seealso \code{\link{create.THMM}}
#'
#' @export
as.TDCM.THMM <- function(x) {
    if (!is.THMM(x)) stop("'x' must be of class 'THMM'")

    # Preparing initial data
    data <- x[c("PTNUM", "time", "state", "covariate")]

    # Initialize a list to store the transformed rows
    transformedRows <- list()

    i <- 1
    while (i <= nrow(data)) {
        currentState <- data$state[i]
        nextState <- ifelse(i < nrow(data), data$state[i + 1], NA)
        nextTime <- ifelse(i < nrow(data), data$time[i + 1], data$time[i])

        # Default auxiliary row setup
        aux <- c(PTNUM = data$PTNUM[i], start = 0, stop = nextTime, event = 0, covariate = data$covariate[i], tdcov = 0)

        # Adjust 'aux' based on state transitions
        if (!is.na(nextState)) {
            if (currentState == 1 && nextState == 3) {
                aux$event <- 1 # Transition from state 1 to 3
            } else if (currentState == 2) {
                aux$tdcov <- 1 # State 2 indicates a time-dependent covariate
                aux$start <- data$time[i - 1] # Adjust start time for tdcov
            }
        }

        # Add to list and move index appropriately
        transformedRows[[length(transformedRows) + 1]] <- aux
        i <- ifelse(currentState == 1 && nextState == 3, i + 2, i + 1)
    }

    # Convert the list to a data frame
    data2 <- do.call(rbind, transformedRows)
    names(data2) <- c("id", "start", "stop", "event", "covariate", "tdcov")
    data2 <- data.frame(data2, stringsAsFactors = FALSE)
    class(data2) <- "TDCM"

    return(data2)
}
