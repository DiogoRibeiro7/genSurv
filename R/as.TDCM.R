is.TDCM <- function(x) {
    is.data.frame(x) && inherits(x, "TDCM")
}


as.TDCM <- function(x) {
    UseMethod("as.TDCM")
}


as.TDCM.default <- function(x) {
    stop("cannot coerce class '", deparse(substitute(x)), "' into class 'TDCM'", domain = NA)
}


as.TDCM.TDCM <- function(x) {
    if (!is.TDCM(x)) stop("'x' must be of class 'TDCM'")
    x
}

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
