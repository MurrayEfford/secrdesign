##############################################################################
## package 'secrdesign'
## GAminnr.R
## adapted from code of Ian Durbach 
## 2022-10-02,06
##############################################################################

# default pen_fn
penfn <- function (traps, sigma) {
    # find out how many detector pairs are between 2.5-3.5 and 3.5-4.5 sigma apart
    breaks <- c(0, 2.499, 3.499, 4.499, Inf) * sigma     # why 4.449? assume typo
    d <- as.matrix(dist(traps))  # for compatibility
    tabulate(cut(d, breaks = breaks))[2:3]
}

#-------------------------------------------------------------------------------
# Objective function
OFenrm <- function (v, 
    alltraps, 
    mask, 
    detectpar,
    detectfn,
    noccasions, 
    detector, 
    D, 
    crit, 
    penalty,
    g_penvector) {
    
    # penalty for too clustered
    if (!is.null(penalty)) {
        penvector <- penalty$pen_fn(alltraps[v,], detectpar$sigma)
        penalty <- penalty$pen_wt * sum(pmax(0, g_penvector - penvector))
    }
    else {
        penalty <- 0
    } 
    
    traps <- subset(alltraps, v)
    
    if (length(detectpar$lambda0) > 1)
        stop ("this implementation does not allow varying lambda0")
    enrm <- Enrm(D = D, traps = traps, mask = mask, noccasions = noccasions, 
        detectpar = detectpar, detectfn = detectfn)
    
    c(-enrm[1], -enrm[2], -enrm[3], penalty-(min(enrm[1],enrm[2])))[crit]    
    
}
#-------------------------------------------------------------------------------

GAminnr <- function(
        mask,
    alltraps, 
    ntraps, 
    detectpar, 
    noccasions,
    detectfn = c("HHN", "HHR", "HEX", "HAN", "HCG"),
    D = NULL,
    penalty = NULL,
    seed = NULL,
    ...){
    
    detectfn <- match.arg(detectfn)
    
    ## criterion to use (1 = En, 2 = Er, 3 = Em, 4 = min(En,Er))
    criterion <- 4   # force min(n,r)
    
    if(missing(mask)) stop("Must supply a 'mask' object (coords of the study area)")
    if(missing(alltraps))   stop("Must supply a 'traps' object (all possible trap locations)")
    
    if (!inherits(mask, "mask")) stop ("mask should be a mask object")
    if (!inherits(alltraps, "traps")) stop ("alltraps should be a traps object")    
    
    detector <- match.arg(detector(alltraps), choices = c("count", "proximity", "multi"))
    if (noccasions == 1 && detector == "multi") stop ("multi detector requires > 1 occasion")
    if(!is.null(seed)) set.seed(seed)
    
    if (ms(mask) || ms(traps)) stop ("mask and traps should be single-session")
    
    #---------------------------------------------------------------------------
    
    if (!is.null(penalty)) {
        # penalty reference vector (Durbach et al. 2021)

        # find distribution of trap spacings on a close to regular grid, to ensure 
        # later optimized grid has spaced enough detectors sufficiently far apart 
        # to get low var(sigma)
        
        # polygon to represent region of interest
        pg <- st_union(gridCells(alltraps))
        # place a grid over the area, with cells pen_gridsigma * sigma apart
        # random origin, detector not material
        cellsize <- penalty$pen_gridsigma * detectpar$sigma
        grid_traps <- make.systematic(region = pg, spacing = cellsize)
        # random starting trap
        xy_rand <- grid_traps[sample.int(nrow(grid_traps), 1), ]  
        dist2pt <- distancetotrap(grid_traps, xy_rand)
        OK <- rank(dist2pt, ties.method = "random") <= ntraps
        # closest ntraps to random start
        grid_traps <- subset(grid_traps, OK)   
        # use default penalty function (see above) if none provided
        if (is.null(penalty$pen_fn)) penalty$pen_fn <- penfn  
        # target vector (e.g., minimum number of traps in each distance bracket)
        g_penvector <- penalty$pen_fn(grid_traps, detectpar$sigma)
    }
    else {
        g_penvector <- NA 
    }
    #---------------------------------------------------------------------------
    
    des <- kofnGA::kofnGA(n = nrow(alltraps), 
        k  = ntraps, 
        OF = OFenrm,
        ...,
        alltraps    = alltraps,
        mask        = mask,
        detectpar   = detectpar,
        noccasions  = noccasions,
        detectfn    = detectfn,
        detector    = detector,
        D           = if (is.null(D)) 1 else D,
        crit        = criterion,
        penalty     = penalty,
        g_penvector = g_penvector
    )
    
    optimaltraps <- subset(alltraps, des$bestsol)
    
    if (!is.null(D)) {
        optimalenrm <- Enrm(D = D, traps = optimaltraps, mask = mask, 
            noccasions = noccasions, detectpar = detectpar, detectfn = detectfn)
    }
    else {
        optimalenrm <- NULL
    }
    
    out <- list(
        mask         = mask, 
        alltraps     = alltraps, 
        detectpar    = detectpar, 
        noccasions   = noccasions,
        detectfn     = detectfn,
        D            = D,
        penalty      = penalty,
        des          = des, 
        optimaltraps = optimaltraps,
        optimalenrm  = optimalenrm
        ## do not include minnrRSE - it depends on extra arguments CF, distribution
    )
    
    class(out) <- "GAminnr"
    out
    
}
