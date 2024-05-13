#' Check if an object is of class CMM
#'
#' This function checks if an object is a data frame and inherits the class "CMM".
#'
#' @param x The object to be checked.
#' @return TRUE if the object is of class CMM, FALSE otherwise.
#' @examples
#' is.CMM(data.frame()) # FALSE
#' is.CMM(CMM()) # TRUE
#' @export
is.CMM <- function(x) {
    is.data.frame(x) && inherits(x, "CMM")
}


#' Convert an object to CMM class
#'
#' This function is a generic method for converting an object to the CMM class.
#' It uses the S3 object-oriented programming system in R.
#'
#' @param x The object to be converted.
#'
#' @return An object of class CMM.
#'
#' @export
as.CMM <- function(x) {
    UseMethod("as.CMM")
}


#' Coerce an object into class 'CMM'
#'
#' This function attempts to coerce an object into class 'CMM'. If the coercion is not possible, an error is thrown.
#'
#' @param x The object to be coerced.
#'
#' @return The coerced object.
#'
#' @examples
#' as.CMM.default(5)
#' # Error: Cannot coerce class '5' into class 'CMM'
#'
#' @export
as.CMM.default <- function(x) {
    stop("Cannot coerce class '", deparse(substitute(x)), "' into class 'CMM'", domain = NA)
}


#' Convert an object to class CMM
#'
#' This function converts an object to class CMM.
#'
#' @param x An object to be converted to class CMM.
#' @return The converted object of class CMM.
#' @export
#'
#' @examples
#' x <- as.CMM.CMM(object)
#' class(x)
#'
#' @seealso
#' \code{\link{is.CMM}}
#'
#' @keywords internal
as.CMM.CMM <- function(x) {
    if (!is.CMM(x)) stop("'x' must be of class 'CMM'")
    x
}


#' Convert TDCM object to CMM object
#'
#' This function converts a TDCM (Time-Dependent Covariate Model) object to a CMM (Continuous-Time Markov Model) object.
#'
#' @param x A TDCM object.
#' @return A CMM object.
#' @export
as.CMM.TDCM <- function(x) {
    if (!is.TDCM(x)) stop("'x' must be of class 'TDCM'")

    # Preparing initial data frame
    data <- x[c("id", "start", "stop", "event", "covariate", "tdcov")]

    # Initialize a list for dynamic row addition
    transformedRows <- list()

    # Iterate over rows of data
    for (i in 1:nrow(data)) {
        currentRow <- data[i, ]
        nextRowExists <- i < nrow(data)
        samePatientNextRow <- nextRowExists && data$id[i] == data$id[i + 1]
        
        # Create transformed rows based on conditions
        baseRow <- list(id = currentRow$id, start = 0, stop = currentRow$stop, 
                        event = currentRow$event, covariate = currentRow$covariate, trans = 1)

        if (currentRow$event == 0) {
            transformedRows[[length(transformedRows) + 1]] <- c(baseRow, trans = 1)
            if (!samePatientNextRow || (samePatientNextRow && data$event[i + 1] == 1)) {
                transformedRows[[length(transformedRows) + 1]] <- c(baseRow, trans = 2)
            }
        } else if (currentRow$event == 1) {
            transformedRows[[length(transformedRows) + 1]] <- c(baseRow, trans = 1)
            if (!samePatientNextRow || (samePatientNextRow && data$event[i + 1] == 0)) {
                transformedRows[[length(transformedRows) + 1]] <- c(baseRow, stop = currentRow$stop, event = 0, trans = 2)
            }
        }

        # Additional case for last row or transitioning patients
        if (!nextRowExists || (nextRowExists && !samePatientNextRow)) {
            transformedRows[[length(transformedRows) + 1]] <- c(baseRow, trans = 3)
        }
    }

    # Convert the list to a data frame
    data2 <- do.call(rbind, transformedRows)
    names(data2) <- c("id", "start", "stop", "event", "covariate", "trans")
    data2 <- data.frame(data2, stringsAsFactors = FALSE)
    class(data2) <- "CMM"
    
    return(data2)
}


#' Convert a THMM object to a CMM object
#'
#' This function converts a THMM (Time-Homogeneous Markov Model) object to a CMM (Continuous-Time Markov Model) object.
#' The THMM object must have the following columns: "PTNUM" (patient ID), "time" (time of observation), "state" (current state), and "covariate" (covariate value).
#' The function creates a new data frame with the following columns: "id" (patient ID), "start" (start time), "stop" (stop time), "event" (event indicator), "covariate" (covariate value), and "trans" (transition indicator).
#' The "event" column indicates whether an event occurred during the time interval, and the "trans" column indicates the type of transition.
#' If a patient has multiple state sequences, the function handles them separately.
#'
#' @param x A THMM object to be converted to a CMM object.
#' @return A CMM object.
#' @export
#' @examples
#' # Create a THMM object
#' thmm <- data.frame(PTNUM = c(1, 1, 1, 2, 2, 2),
#'                    time = c(0, 1, 2, 0, 1, 2),
#'                    state = c(1, 2, 3, 1, 2, 3),
#'                    covariate = c(0.5, 0.7, 0.9, 1.2, 1.5, 1.8))
#' class(thmm) <- "THMM"
#'
#' # Convert the THMM object to a CMM object
#' cmm <- as.CMM.THMM(thmm)
#' class(cmm)
#' # Output: "CMM"
as.CMM.THMM <- function(x) {
    if (!is.THMM(x)) stop("'x' must be of class 'THMM'")
    
    # Preparing the initial data frame directly
    data <- x[c("PTNUM", "time", "state", "covariate")]

    # Initialize a list to collect rows for efficiency
    transformedRows <- list()
    
    i <- 1
    while (i <= nrow(data)) {
        patientID <- data$PTNUM[i]
        covariate <- data$covariate[i]
        
        # Check for start of a new patient's data or a new state sequence
        if (data$time[i] == 0) {
            nextState <- ifelse(i + 1 <= nrow(data) && data$PTNUM[i + 1] == patientID, data$state[i + 1], NA)
            nextTime <- ifelse(!is.na(nextState), data$time[i + 1], data$time[i])
            
            # Creating auxiliary rows based on the next state
            auxRows <- list(
                c(patientID, 0, nextTime, nextState %in% c(NA, 3), covariate, 1),
                c(patientID, 0, nextTime, nextState == NA || nextState == 3, covariate, 2)
            )
            
            # Handle specific cases based on nextState and potentially skip rows
            if (!is.na(nextState) && nextState == 2 && (i + 2 <= nrow(data))) {
                # If there's a following state 2 or 3, create an additional transition
                thirdState <- data$state[i + 2]
                thirdTime <- data$time[i + 2]
                auxRows[[3]] <- c(patientID, nextTime, thirdTime, thirdState == 3, covariate, 3)
                i <- i + 2  # Skip ahead as we've handled two transitions
            }
            
            # Add generated rows to the list
            transformedRows <- c(transformedRows, auxRows)
        }
        
        i <- i + 1
    }
    
    # Convert the list to a data frame
    data2 <- do.call(rbind, transformedRows)
    colnames(data2) <- c("id", "start", "stop", "event", "covariate", "trans")
    data2 <- data.frame(data2, stringsAsFactors = FALSE)
    class(data2) <- "CMM"
    
    return(data2)
}

