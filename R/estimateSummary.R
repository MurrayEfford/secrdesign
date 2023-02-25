
# Convert output (nested list of estimate tables) to array
# with dimensions (Parameter, statistic, Group, Scenario, Replicate)
estimateArray <- function (object) {
    nscen <- length(object$output)
    nrepl <- length(object$output[[1]])
    groups <- unique(object$scenarios$group)
    ngrp  <- max(1, length(groups))
    
    typical <- object$output[[1]][[1]]
    if(!is.data.frame(typical)) {
        groupnames <- if (length(typical) == ngrp) groups else names(typical)
        ngrp <- length(groupnames)
        typical <- typical[[1]]
    }
    else {
        ngrp <- 1
        groupnames <- 1
    }
    num <- sapply(typical, is.numeric)
    if (inherits(object, 'estimatetables')) {
        values <- unlist(object$output)
        values <- suppressWarnings(as.numeric(values))  # link text -> NA
        arr <- array(values, 
            dim = c(dim(typical), ngrp, nrepl, nscen),
            dimnames = list(
                Parameter = rownames(typical),
                names(typical),
                Group = groupnames,
                Replicate = 1:nrepl,
                Scenario = 1:nscen))
        arr <- arr[,num,,, ,drop = FALSE]   # drop link field
        arr <- aperm(arr, c(1,2,3,5,4))
        arr
        
    }
    else stop("unsuitable data; expecting estimatetables from secrdesign")
}

# Direct summary of output from run.scenarios
estimateSummary <- function (object, parameter = 'D', 
    statistics = c('n', 'estimate', 'SE.estimate', 'RB', 'seRB', 
        'RSE', 'RMSE', 'COV'), true, validrange = c(0, Inf), 
    checkfields = c('estimate','SE.estimate'),
    format = c('list', 'data.frame')) {
    
    format <- match.arg(format)
    
    if (missing(true)) {
        true <- object$scenarios[,parameter]
    }
    arr <- estimateArray(object)

    if (length(true) != prod(dim(arr)[3:4])) stop ("incongruent 'true'")
    # check for out-of-range estimates
    validate <- function(x) {
        x1 <- as.numeric(x[checkfields])
        if (!all(x1 >= validrange[1] & x1 <= validrange[2])) {
            x[] <- NA
        }
        x
    }
    arr <- aperm(apply(arr, c(1,3,4,5), validate), c(2,1,3,4,5))
    
    # subsets of data
    parm <- arr[parameter,'estimate',,,, drop = FALSE]
    parmse <- arr[parameter,c('estimate','SE.estimate'),,,, drop = FALSE]
    parmcl <- arr[parameter, c('lcl','ucl'),,,, drop = FALSE]
    
    se <- function (x) sd(x, na.rm = TRUE) / sum(!is.na(x))^0.5
    n <- apply(parm, 3:4, function(x) sum(!is.na(x)))
    EST <- apply(parm, 3:4, mean, na.rm = TRUE)
    seEST <- apply(parm, 3:4, se)
    
    # relative bias
    rb <- function (x, true) (x-true)/true
    RB <- sweep(parm, MARGIN = 3:4, STATS = true, FUN = rb)
    mnRB <- apply(RB, 3:4, mean, na.rm = TRUE)
    seRB <- apply(RB, 3:4, se)
    
    # relative SE
    rse <- function (x) x[2]/x[1]
    RSE <- apply(parmse, 3:5, rse)
    RSE <- apply(RSE, 1:2, mean, na.rm = TRUE)
    
    # relative RMSE
    err <- function (x, true) (x-true)
    rms <- function (x) sqrt(mean(x^2, na.rm = TRUE))
    ERR <- sweep(parm, MARGIN = 3:4, STATS = true, FUN = err)
    RMSE <- apply(ERR, 3:4, rms) / true
    
    # coverage of CI
    inrange <- function(x) (x[1,1,,,drop=FALSE] < true) & (x[1,2,,,drop=FALSE] > true)
    OK <- apply(parmcl, 5, inrange)
    dim(OK) <- dim(parmcl)[3:5]
    COV <- apply(OK, 1, mean, na.rm = TRUE)
    COV <- array(COV, dim = dim(parm)[3:4], dimnames = dimnames(parm)[3:4])
    
    out <- list(n = n, estimate = EST, SE.estimate = seEST, RB = mnRB, 
        seRB = seRB, RSE = RSE, RMSE = RMSE, COV = COV)[statistics]
    if (format == 'list') 
        out
    else
        do.call(data.frame, lapply(out, as.numeric))
}
