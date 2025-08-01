##############################################################################
## package 'secrdesign'
## make.scenarios.R
## 2014-02-06, 2014-02-09

## 2014-02-19 popindex and popargs replace Ndist; FIXED 'CROSSING'
## 2014-04-27 detindex
## 2014-11-25 groups

###############################################################################

## construct a dataframe in which each row represents a scenario
make.scenarios <- function (trapsindex = 1, noccasions = 3, nrepeats = 1,
                            D, 
                            g0, sigma, lambda0, epsilon, tau, z,
                            detectfn = 0, recapfactor = 1,
                            popindex = 1, detindex = 1, fitindex = 1, groups,
                            crosstraps = TRUE) {
    inputs <-  as.list (environment())
    inputs$crosstraps <- NULL
    if (missing(D)) {
        inputs$nrepeats <- NA
        inputs$D <- NA
    }
    trapping   <- inputs[c('trapsindex', 'noccasions', 'nrepeats')]
    
    ## allow uniform detectfn = 4 and detectfn = 20 for simulation
    inputs$detectfn <- secr:::secr_valid.detectfn(inputs$detectfn, valid = c(0:20))
    pnames <- secr:::secr_parnames(inputs$detectfn)
    OK <- sapply(inputs[pnames], is.numeric)
    if (any (!OK)) stop ("parameters missing from input: ", pnames[!OK])
    
    parameters <- inputs[c('D', pnames, 'detectfn', 'recapfactor',
                           'popindex', 'detindex', 'fitindex')]
    
    if (!crosstraps) {
        trapmat <- matrix(nrow = max(sapply(trapping, length)), ncol = 3)
        for (i in 1:3) trapmat[,i] <- trapping[[i]]
        trapdf <- data.frame(trapmat)
        names(trapdf) <- c('trapsindex', 'noccasions', 'nrepeats')
        trapping <- list(trapping = 1:nrow(trapmat))
    }
    if (!all(sapply(parameters[-1], is.numeric)))
        ## stop ("must provide all detection parameters")
        warning ("not all detection parameters provided: complete manually")

    value <- do.call (expand.grid, c(trapping, parameters, list(stringsAsFactors = FALSE)))

    ## repeat to fill column if trap settings are not crossed
    if (!crosstraps) {
        ## assume trapping index in col 1
        value <- cbind(trapdf[value[,1],], value[,-1])
    }
    nr <- nrow(value)
    value <- cbind (scenario = 1:nr, value)
    if (!missing(groups)) {
        # groups <- factor(groups)
        ng <- length(groups)
        value <- cbind(value[rep(1:nr, each = ng),], group = rep(groups, nr))[,c(1,13,2:12)]
    }
    rownames(value) <- 1:nrow(value)
    attr(value, 'inputs') <- inputs    ## used in make.array()
    value
}

# make.scenarios(D = 1:2, g0 = seq(0.1,0.3,0.1), sigma=25)
# make.scenarios(D = 1:2, sigma=25, lambda0 = 1)
# make.scenarios(trap=1:3, nrepeats=3:1, D = 1:4, sigma=25, lambda0 = 1, cross=F)
# tmp <- attr(make.scenarios(D = 1:2, g0 = seq(0.1,0.3,0.1), sigma = 25), 'inputs')
# do.call(make.scenarios, tmp)

