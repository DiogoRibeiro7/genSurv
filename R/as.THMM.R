#' Check if an object is of class THMM
#'
#' This function checks if an object is of class THMM by verifying if it is a data frame and if it inherits the "THMM" class.
#'
#' @param x The object to be checked
#' @return TRUE if the object is of class THMM, FALSE otherwise
#' @examples
#' is.THMM(data.frame()) # FALSE
#' is.THMM(THMM()) # TRUE
is.THMM <- function(x) {
    is.data.frame(x) && inherits(x, "THMM")
}


#' Convert an object to THMM class
#'
#' This function is a generic method for converting an object to the THMM class.
#' It dispatches the appropriate method based on the class of the object.
#'
#' @param x The object to be converted.
#' @return An object of class THMM.
#' @export
as.THMM <- function(x) {
    UseMethod("as.THMM")
}


#' Coerce an object into class 'THMM'
#'
#' This function is used to coerce an object into class 'THMM'. If the coercion is not possible, an error is thrown.
#'
#' @param x The object to be coerced.
#' @return An object of class 'THMM'.
#' @export
as.THMM.default <- function(x) {
    stop("cannot coerce class '", deparse(substitute(x)), "' into class 'THMM'", domain = NA)
}


#' Convert object to THMM class
#'
#' This function converts an object to the THMM class.
#'
#' @param x The object to be converted.
#' @return The object of class THMM.
#' @export
#' @examples
#' as.THMM.THMM(object)
as.THMM.THMM <- function(x) {
    if (!is.THMM(x)) {
        stop("'x' must be of class 'THMM'")
    }
    x
}


#' Convert a CMM object to a THMM object
#'
#' This function converts a CMM (Continuous-time Markov Model) object to a THMM (Two-state Hidden Markov Model) object.
#' The CMM object should have the following columns: "id", "start", "stop", "event", "covariate", and "trans".
#' The THMM object will have the following columns: "PTNUM", "time", "state", and "covariate".
#'
#' @param x A CMM object to be converted to a THMM object.
#' @return A THMM object.
#' @export
as.THMM.CMM <- function(x) {
    if (!is.CMM(x)) stop("'x' must be of class 'CMM'")

    # Simplify the data frame creation with direct assignment
    data <- x[, c("id", "start", "stop", "event", "covariate", "trans")]
    names(data) <- c("id", "start", "stop", "event", "covariate", "trans")

    # Initialize an empty list for storing transformed rows
    transformedRows <- list()

    i <- 1
    while (i <= nrow(data) - 1) {
        currentState <- data[i, "event"]
        nextState <- data[i + 1, "event"]
        patientID <- data[i, "id"]
        covariate <- data[i, "covariate"]

        if (currentState == 0 && nextState == 0) {
            transformedRows[[length(transformedRows) + 1]] <- c(patientID, 0, 1, covariate)
            transformedRows[[length(transformedRows) + 1]] <- c(patientID, data[i, "stop"], 1, covariate)
            i <- i + 2
        } else if (currentState == 1 && data[i, "trans"] == 1) {
            transformedRows[[length(transformedRows) + 1]] <- c(patientID, 0, 1, covariate)
            transformedRows[[length(transformedRows) + 1]] <- c(patientID, data[i, "stop"], 3, covariate)
            i <- i + 2
        } else if (currentState == 0 && nextState == 1) {
            nextStateDuration <- data[min(i + 2, nrow(data)), "stop"]
            nextStateEvent <- data[min(i + 2, nrow(data)), "event"]

            transformedRows[[length(transformedRows) + 1]] <- c(patientID, 0, 1, covariate)
            transformedRows[[length(transformedRows) + 1]] <- c(patientID, data[i, "stop"], 2, covariate)
            transformedRows[[length(transformedRows) + 1]] <- c(patientID, nextStateDuration, nextStateEvent == 0 ? 2:3, covariate)

            i <- i + 3 # Assume the next state is always present
        }
    }

    # Convert the list of rows into a data frame
    data2 <- do.call(rbind, transformedRows)
    colnames(data2) <- c("PTNUM", "time", "state", "covariate")

    data2 <- data.frame(data2, stringsAsFactors = FALSE, row.names = NULL)
    row.names(data2) <- as.integer(1:nrow(data2))
    class(data2) <- c("data.frame", "THMM")

    return(data2)
}


#' Convert TDCM object to THMM object
#'
#' This function converts a TDCM (Time-Dependent Covariate Model) object to a THMM (Time-Homogeneous Markov Model) object.
#'
#' @param x A TDCM object.
#' @return A THMM object.
#' @export
as.THMM.TDCM <- function(x) {
    if (!is.TDCM(x)) stop("'x' must be of class 'TDCM'")

    # Initial setup
    data <- x[c("start", "stop", "event", "covariate", "tdcov")]
    transformedData <- list()
    patientIndex <- 1

    # Loop through data to transform
    for (i in 1:nrow(data)) {
        currentState <- data$event[i]
        nextRowExists <- i < nrow(data)
        transitionOccursNext <- nextRowExists && (data$start[i + 1] == data$stop[i])
        eventIsCensoredNext <- transitionOccursNext && (data$event[i + 1] == 0)

        # Common row to add for each patient
        commonRow <- list(PTNUM = patientIndex, time = 0, state = 1, covariate = data$covariate[i])

        if (eventIsCensoredNext) {
            # Handling continuous censored transitions
            transformedData[[length(transformedData) + 1]] <- commonRow
            transformedData[[length(transformedData) + 1]] <- c(commonRow$PTNUM, data$stop[i], 2, commonRow$covariate)
            i <- i + 1 # Skip the next row as it has been handled
        } else {
            nextState <- ifelse(currentState == 1, 3, 1)
            if (nextRowExists) nextState <- ifelse(data$event[i + 1] == 1, 3, 2)

            transformedData[[length(transformedData) + 1]] <- commonRow
            if (nextState != 1 || !eventIsCensoredNext) {
                stopTime <- ifelse(nextRowExists, data$stop[i + 1], data$stop[i])
                transformedData[[length(transformedData) + 1]] <- c(commonRow$PTNUM, data$stop[i], nextState, commonRow$covariate)
            }

            if (nextState == 3 && nextRowExists) i <- i + 1 # Additional skip for directly observed event
        }

        patientIndex <- patientIndex + 1
        if (!nextRowExists) break # Exit loop if at the end
    }

    # Finalize the data frame
    data2 <- do.call(rbind, transformedData)
    names(data2) <- c("PTNUM", "time", "state", "covariate")
    data2 <- data.frame(data2, row.names = NULL, stringsAsFactors = FALSE)
    class(data2) <- c("data.frame", "THMM")

    return(data2)
}
