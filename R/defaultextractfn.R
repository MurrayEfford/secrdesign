defaultextractfn <- function(x, ...) {
    
    if (inherits(x, 'try-error')) {
        ## null output: dataframe of 0 rows and 0 columns
        data.frame()
    }
    else if (inherits(x, 'capthist')) {
        ## summarised raw data
        if (ms(x))
            unlist(lapply(x, counts))
        else {
            gp <- covariates(x)$group
            if (is.null(gp))
                counts(x)
            else
                unlist(lapply(split(x,gp,dropnullocc=TRUE), counts))
        }
    }
    else if (
        (inherits(x, 'secr') && !is.null(x$fit)) ||
        (inherits(x, 'ipsecr') && x$code == 1)
    ) {
        ## fitted model:
        ## default predictions of 'real' parameters
        out <- predict(x)
        if (!is.data.frame(out)) {
            warning ("summarising only first session, group or mixture class")
            out <- out[[1]]
        }
        attr(out, 'counts')   <- counts(x$capthist)
        attr(out, 'start')    <- x$start
        attr(out, 'method')   <- x$method
        attr(out, 'proctime') <- x$proctime
        attr(out, 'chat')     <- x$details$chat
        attr(out, 'D')        <- attr(x$capthist, 'D')    # usually NULL
        
        # ipsecr.fit completion 'code' 1 successful, 2 target not within final box, 3 exceeded maximum simulations
        # nlm() completion 'code' 
        # optim() completion code 'convergence'
        if (!is.null(x$fit) && !is.null(x$fit$code))
            attr(out, 'code')   <- x$fit$code
        if (!is.null(x$fit) && !is.null(x$fit$convergence))
            attr(out, 'convergence') <- x$fit$convergence
        
        out
    }
    else {
        ## null output: dataframe of 0 rows and 0 columns
        data.frame()
    }
}

