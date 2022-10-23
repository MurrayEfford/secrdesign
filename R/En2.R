##############################################################################
## package 'secrdesign'
## En2.R
## 2022-10-23
##############################################################################

En2 <- function (D, traps, mask, detectpar, noccasions, detectfn = 
        c('HHN', 'HHR', 'HEX','HAN','HCG', 'HN', 'HR', 'EX')) {
    
    if (is.character(detectfn))
        detectfn <- match.arg(detectfn)
    detectfn <- secr:::valid.detectfn(detectfn, valid = c(0,1,2,14:18))
    dfc <- dfcast (detectfn, detectpar)  # transforms detectfn 0 to 14, 2 to 16
    detectfn <- dfc$detectfn
    detectpar <- dfc$detectpar
    detectpars <- unlist(detectpar[secr:::parnames(detectfn)])
    dettype <- secr:::detectorcode(traps, noccasions = noccasions)[1]
    D <- rep(D, length.out = nrow(mask)) * attr(mask, 'area')  # per cell
    temp <- En2cpp (
        as.integer(dettype),
        unlist(detectpars), 
        as.double(D),
        edist(traps, mask), 
        as.integer(detectfn),
        as.integer(noccasions)
    )
    
    c(En = temp$En, En2 = temp$En2)
    
}