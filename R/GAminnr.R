# adapted from code of Ian Durbach 
# 2022-10-02,4

# replace 'statespace' with 'mask'
# replace e2dist with secr function edist
# remove commented code (#N etc.)
# replace "count" by detector in read.traps()
# in Appendix F D_per_mask_cell is not an argument, use D = D,
# share function td for penalty term
# return more variables

# all_grid_traps <- list() removed

# remove unused arguments from OFenrm: 
# pop, N = 100, sum_wt, transitions = NULL, divby = 1, addback = 0 

# USE THE FOLLOWING WITH CAUTION!
# pen_gridsigma:   penalties for bad trap spacing are based on comparing min(n,r) grid
#                  to what you would get from a X-sigma regular grid. X = pen_gridsigma
# pen_wt     :          penalty for spacings that are too different from regular grid

#-------------------------------------------------------------------------------

td <- function (traps, sigma) {
    # find out how many detector pairs are between 2-3 and 3-4 sigma apart
    breaks <- c(0, 1.499, 2.499, 3.499, 4.449, Inf) * sigma
    d <- as.matrix(dist(traps))  # for compatibility
    tabulate(cut(d, breaks = breaks))[2:4]
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
    g_n234, 
    pen_wt,
    pen_fn) {
    
    lambda0 <- detectpar$lambda0
    sigma <- detectpar$sigma
    
    # penalty for too clustered
    if (pen_wt > 0) {
        # use function passed to OF
        n234 <- pen_fn(alltraps[v,], detectpar$sigma)
        penalty <- pen_wt * sum(pmax(0, g_n234-n234))
    }
    else {
        penalty <- 0
    } 
    
    traps <- subset(alltraps, v)
    
    # this implementation does not allow varying lambda0
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
    noccasions = 1,
    detectfn = c("HHN", "HHR", "HEX", "HAN", "HCG"),
    D = 1,
    criterion = 4,     # default min(n,r)
    pen_gridsigma = 2,
    pen_wt = 0, 
    seed = NULL,
    ...){
    
    detectfn <- match.arg(detectfn)
    
    if(missing(mask)) stop("Must supply a 'mask' object (coords of the study area)")
    if(missing(alltraps))   stop("Must supply a 'traps' object (all possible trap locations)")
    
    if (!inherits(mask, "mask")) stop ("mask should be a mask object")
    if (!inherits(alltraps, "traps")) stop ("alltraps should be a traps object")    
    
    detector <- match.arg(detector(alltraps), choices = c("count", "proximity", "multi"))
    if(!is.null(seed)) set.seed(seed)
    
    if (ms(mask) || ms(traps)) stop ("mask and traps should be single-session")
    
    #---------------------------------------------------------------------------
    # penalty base
    
    if (pen_wt>0) {
        # find distribution of trap spacings on a close to regular grid, to ensure later optimized grid has spaced
        # enough detectors sufficiently far apart to get low var(sigma)
        
        pg <- st_union(gridCells(alltraps))
        # place a grid over the area, with cells X sigma apart
        cellsize <- pen_gridsigma * detectpar$sigma
        grid_traps <- make.systematic(region = pg, spacing = cellsize, detector = detector)  # random origin
        xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
        dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
        OK <- rank(dist2pt, ties.method = "random") <= ntraps
        grid_traps <- subset(grid_traps, OK)
        
        # target for minimum numbers of traps in each distance bracket
        targetmultiplier <- c(0,1,1)   # 0.8, 1, 1 ?
        g_n234 <- td(grid_traps, detectpar$sigma) * targetmultiplier
    }
    else {
        g_n234 <- c(1,1,1) 
    }
    #---------------------------------------------------------------------------
    
    des <- kofnGA::kofnGA(n = nrow(alltraps), 
        k  = ntraps, 
        OF = OFenrm,
        ...,
        alltraps   = alltraps,
        mask       = mask,
        detectpar  = detectpar,
        noccasions = noccasions,
        detectfn   = detectfn,
        detector   = detector,
        D          = D,
        crit       = criterion,
        pen_wt     = pen_wt,
        pen_fn     = td,
        g_n234     = g_n234
    )
    
    optimaltraps <- subset(alltraps, des$bestsol)
    
    optimalenrm <- Enrm(D = D, traps = optimaltraps, mask = mask, 
        noccasions = noccasions, detectpar = detectpar, detectfn = detectfn)
    
    list(
        mask         = mask, 
        alltraps     = alltraps, 
        detectpar    = detectpar, 
        noccasions   = noccasions,
        detectfn     = detectfn,
        D            = D,
        des          = des, 
        optimaltraps = optimaltraps,
        optimalenrm  = optimalenrm
    )
    
}
