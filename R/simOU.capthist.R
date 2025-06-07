simOU <- function(xy, tau, sigma, noccasions, start = NULL){
    xy <- as.numeric(xy)   # dodge issue with dataframe
    if (!requireNamespace("mvtnorm")) stop("simOU requires package mvtnorm")
    if (is.null(start)) {
        start <- mvtnorm::rmvnorm(1, mean = xy, sigma = sigma^2*diag(2))
    }
    beta <- 1/tau
    out <- matrix(0, nrow = noccasions, ncol = 2)
    out[1, ] <- start
    for (i in 2:noccasions){
        out[i, ] <- mvtnorm::rmvnorm(1, 
            mean = (1 - exp(-beta)) * xy + exp(-beta) * out[i - 1, ],
            sigma = sigma^2 * (1 - exp(-2*beta)) * diag(2))
    }
    out
}

simOU.capthist <- function (
        traps,
        popn,
        detectpar,     # list of tau/lambda0, sigma
        noccasions,    # effective "duration"
        epsilon,       # proximity threshold radius for detection
        savepopn = FALSE,
        savepath = FALSE,
        ...)
{
    captfn <- function (xy) secr::edist(xy,traps) <= epsilon   
    N <- nrow(popn)
    tau <- detectpar$tau
    if (is.null(tau)) tau <- detectpar$lambda0
    locs <- apply(popn, 1, simOU, tau, detectpar$sigma, noccasions, simplify = FALSE)
    capt <- lapply(locs, captfn)
    capt <- do.call(rbind, capt)  
    ch <- array(capt, dim = c(noccasions, N, ncol(capt)), 
                dimnames = list(1:noccasions,rownames(popn),rownames(traps)))
    ch <- aperm(ch, c(2,1,3))
    ch <- ch[apply(ch, 1, sum) > 0,,, drop = FALSE]   # drop null histories
    class(ch) <- 'capthist'
    traps(ch) <- traps
    # cast as required detector type
    ch <- reduce(ch, outputdetector = detector(traps)[1], dropunused = FALSE, ...)
    if (savepopn) attr(ch, 'popn') <- popn
    if (savepath) attr(ch, 'path') <- locs
    ch
}
