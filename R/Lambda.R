##############################################################################
## package 'secrdesign'
## Lambda.R
## 2017-05-04
## 2017-07-15 getRSE renamed minnrRSE; optimalSpacing moved to optimalSpacing.R
## 2017-07-21 getdetectpar moved to getdetectpar.R
##############################################################################

Lambda <- function (traps, mask, detectpar, noccasions, detectfn = 
                        c('HHN', 'HHR', 'HEX','HAN','HCG', 'HN', 'HR', 'EX')) {
    if (ms(traps) | ms(mask))
        stop ("Lambda does not accept multi-session traps or mask")
    if (is.character(detectfn))
        detectfn <- match.arg(detectfn)
    detectfn <- secr:::valid.detectfn(detectfn, valid = c(0,1,2,14:18))
    truncate <- ifelse(is.null(detectpar$truncate), 1e+10, detectpar$truncate)
    dfc <- dfcast (detectfn, detectpar)  # transforms detectfn 0 to 14, 2 to 16
    detectfn <- dfc$detectfn
    detectpar <- dfc$detectpar
    detectpars <- unlist(detectpar[secr:::parnames(detectfn)])
    dettype <- secr:::detectorcode(traps, noccasions = noccasions)
    
    temp <- .C("LambdaC", 
        as.double(detectpars),
        as.integer(nrow(traps)),
        as.integer(nrow(mask)),
        as.double(unlist(traps)),
        as.double(unlist(mask)),
        as.integer(detectfn),
        L = double(nrow(mask)),
        resultcode = integer(1))
    
    if (temp$resultcode != 0)
        stop ("error in external function 'LambdaC'")
    L <- temp$L * noccasions
    covariates(mask) <- data.frame(Lambda = temp$L * noccasions)
    
    temp <- .C("sumpkC", 
               as.integer(dettype),
               as.double(detectpars),
               as.integer(nrow(traps)),
               as.integer(nrow(mask)),
               as.double(unlist(traps)),
               as.double(unlist(mask)),
               as.integer(detectfn),
               sumpk = double(nrow(mask)),
               sumq2 = double(nrow(mask)),
               resultcode = integer(1))
    if (temp$resultcode != 0)
        stop ("error in external function 'sumpkc'")
    covariates(mask)$sumpk <- temp$sumpk * noccasions
    covariates(mask)$sumq2 <- temp$sumq2 
    return(mask)
}

Lambdak <- function (D, traps, mask, detectpar, detectfn = 
        c('HHN', 'HHR', 'HEX','HAN','HCG', 'HN', 'HR', 'EX')) {
    if (ms(traps) | ms(mask))
        stop ("Lambdak does not accept multi-session traps or mask")
    if (is.character(detectfn))
        detectfn <- match.arg(detectfn)
    detectfn <- secr:::valid.detectfn(detectfn, valid = c(0,1,2,14:18))
    truncate <- ifelse(is.null(detectpar$truncate), 1e+10, detectpar$truncate)
    dfc <- dfcast (detectfn, detectpar)  # transforms detectfn 0 to 14, 2 to 16
    detectfn <- dfc$detectfn
    detectpar <- dfc$detectpar
    detectpars <- unlist(detectpar[secr:::parnames(detectfn)])
    if (any(detector(traps) != 'capped')) warning ("Lambdak is intended only for capped detectors")
    dettype <- 8
    temp <- .C("LambdaK", 
        as.double(detectpars),
        as.integer(nrow(traps)),
        as.integer(nrow(mask)),
        as.double(unlist(traps)),
        as.double(unlist(mask)),
        as.integer(detectfn),
        LK = double(nrow(traps)),
        resultcode = integer(1))
    
    if (temp$resultcode != 0)
        stop ("error in external function 'LambdaK'")
    return( D * temp$LK * attr(mask,'area'))  
}

Encap <- function (D, traps, mask, detectpar, noccasions, detectfn, LK) {
    detectfn <- secr:::valid.detectfn(detectfn, valid = c(0,1,2,14:18))
    detectpars <- unlist(detectpar[secr:::parnames(detectfn)])
    if (any(detector(traps) != 'capped')) warning ("Lambdak is intended only for capped detectors")
    dettype <- 8

    temp <- .C("Encapped", 
        as.double(detectpars),
        as.integer(nrow(traps)),
        as.integer(nrow(mask)),
        as.double(unlist(traps)),
        as.double(unlist(mask)),
        as.integer(detectfn),
        as.integer(noccasions),
        as.double(LK),
        as.double(D * attr(mask, 'area')),
        En = double(1),
        resultcode = integer(1))
    
    if (temp$resultcode != 0)
        stop ("error in external function 'LambdaK'")
    return(temp$En)  
}

Enrm <- function (D, ...) {
    allargs <- list(...)
    tr <- which(sapply(allargs, inherits, "traps"))
    detect <- detector(allargs[[tr]])[1]
    if (!detect %in% c("multi","proximity","count","capped"))
        stop("only for 'multi', 'proximity', 'count' and 'capped' detectors")
    cellarea <- attr(allargs$mask, 'area')

    if (detect == "capped") {
        warning ("expected values are biased for capped detectors because of spatial heterogeneity")
        nocc <- match("noccasions", names(allargs))
        LK <- do.call(Lambdak, c(list(D), allargs[-nocc]))
        Ec <- allargs$noccasions * sum(1 - exp(-LK))
        En <- do.call(Encap, c(list(D = D), allargs, list(LK = LK)))
        Em <- NA
    }
    else {
        L <- Lambda(...)
        Lam <- covariates(L)$Lambda
        En <- sum( D * cellarea * (1 - exp(-Lam))) # expected n; multi, proximity, count
        if (detect %in% "count") { 
            Ec <- sum(D * cellarea * Lam )
            Em <- sum(D * cellarea * (covariates(L)$Lambda - (1-exp(-covariates(L)$Lambda))) * 
                    (1-covariates(L)$sumq2))
        }
        else {   ## detect %in% c("multi", "proximity")
            Ec <- sum(D * cellarea * covariates(L)$sumpk )
            Em <- sum(D * cellarea * ((covariates(L)$sumpk - (1-exp(-covariates(L)$Lambda))) * 
                    (1-covariates(L)$sumq2)))
        }
    }
    Er <- Ec - En        # expected r; count
    c(En = En, Er = Er, Em = Em)
}

minnrRSE <- function (D, ..., CF = 1.0, distribution = c('poisson', 'binomial')) {
    distribution <- match.arg(distribution)
    nrm <- Enrm(D, ...)
    RSE <- sqrt(CF/min(nrm[1:2]))
    if (distribution == 'binomial') {
        mask <- list(...)$mask
        # if not named, assume second
        if (is.null(mask)) mask <- list(...)[[2]]
       A <- maskarea(mask) 
       RSE <- sqrt (RSE^2 - 1 / (D * A))
    }
    RSE
}