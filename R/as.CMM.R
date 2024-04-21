is.CMM <- function(x) {
    is.data.frame(x) && inherits(x, "CMM")
}


as.CMM <- function(x) {
    UseMethod("as.CMM")
}


as.CMM.default <- function(x) {
    stop("Cannot coerce class '", deparse(substitute(x)), "' into class 'CMM'", domain = NA)
}


as.CMM.CMM <- function(x) {
    if (!is.CMM(x)) stop("'x' must be of class 'CMM'")
    x
}


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

