###############################################################################
## package 'secrdesign'
## select.stats.R
## 2014-11-25 moved from methods.R
## 2015-02-17 user-specified 'true' values
## 2016-02-04 fix weighted problem for param != 'D'
###############################################################################

weighted <- function (onescenario, param) {
  ## 2016-02-04 rename argument scenario to onescenario to ensure OK for param != 'D'
    with(onescenario, {
        if (param == 'D')
            sum(D)
        else if (param == 'pmix') {
            warning("assuming first scenario row corresponds to requested pmix")
            D[1]/sum(D)
        }
        else {
            wt <- D/sum(D)
            sum(onescenario[,param] * wt)
        }
    })
}

select.stats <- function (object, parameter = 'D', statistics, true) {
  
  if (!inherits(object, 'estimatetables'))
    stop ("select.stats requires input of class estimatetables")
  if (is.na(object$outputtype))
    stop ("cannot select.stats output of unknown type")
  if (object$outputtype %in% c('secrfit'))
    stop ("cannot select.stats from fitted model - use predict() first")
  #if (missing(estname) | missing(SEname)) {
  estname <- ""
  SEname <- ""
  if (object$outputtype %in% c('predicted', 'derived', 'regionN')) {
      estname <- 'estimate'
      SEname <- 'SE.estimate'
  }
  else if (object$outputtype %in% c('rawcounts')) {
        estname <- 'n'
        SEname <- 'r'
        parameter <- 'Number'
    }
    else if (object$outputtype == 'coef') {
        estname <- 'beta'
        SEname <- 'SE.beta'
    }
    #}
    for (i in 1:length(object$output)) {
        # typical <- object$output[[i]][[1]]  ## ith scenario, first replicate
        typical <- data.frame(object$output[[i]][[1]])  ## ith scenario, first replicate
        if (length(typical) > 0) break
    }
    if (length(typical) == 0) stop ("no results found")
 
    stat0 <- names(typical)[sapply(typical, is.numeric)]
    if (missing(statistics)) {
        stat1 <- stat0
        if (parameter == 'Number') 
            stat2 <- character(0)
        else
            stat2 <- c('RB','RSE','COV')
        ## stat2 <- c('true','RB','RSE','COV','ERR')
    }
    else {
        stat1 <- statistics[statistics %in% stat0]
        stat2 <- statistics[statistics %in% c('true','RB','RSE','COV','ERR')]
    }
    if (any(stat2 %in% c('true','RB','RSE','COV','ERR')) & (estname == '')) {
        stat2 <- character(0)
        warning ("cannot compute all requested statistics with your data")
    }

    extractfn <- function (out, true, estimated) {
        getij <- function(df, i, j) {
            if (nrow(df) == 0)
                rep(NA, length(j))
            else
                df[i,j]
        }
        if (length(stat1)>0) {
            tmp <- lapply(out, getij, parameter, stat1)
            tmp <- do.call (rbind, tmp)
            rownames(tmp) <- 1:nrow(tmp)
        }
        else
            tmp <- matrix(nrow = length(out), ncol = 0)

        for (st in stat2) {
            if (st == 'true') {
                tmp <- cbind(tmp, rep(true, nrow(tmp)))
            }
            if (st == 'RB') {
                if (estimated) {
                    est <- sapply(out, getij, parameter, estname)
                    tmp <- cbind(tmp, (est - true) / true)
                }
                else tmp <- cbind (tmp, rep(NA, nrow(tmp)))
            }
            if (st == 'RSE') {
                est <- sapply(out, getij, parameter, estname)
                SE.est <- sapply(out, getij, parameter, SEname)
                tmp <- cbind(tmp, SE.est/est)
            }
            if (st == 'ERR') {
                if (estimated) {
                    est <- sapply(out, getij, parameter, estname)
                    tmp <- cbind(tmp, abs(est-true))
                }
                else tmp <- cbind (tmp, rep(NA, nrow(tmp)))
            }
            if (st == 'COV') {
                if (estimated) {
                    est <- sapply(out, getij, parameter, estname)
                    lcl <- sapply(out, getij, parameter, 'lcl')
                    ucl <- sapply(out, getij, parameter, 'ucl')
                    tmp <- cbind(tmp, as.numeric((true>lcl) & (true<ucl)))
                }
                else tmp <- cbind (tmp, rep(NA, nrow(tmp)))
            }
        }
        colnames(tmp) <- c(stat1, stat2)
        tmp
    }

    uniqueScenarioIndex <- match(unique(object$scenarios$scenario), object$scenarios$scenario)
    fitIndex <- object$scenarios$fitindex[uniqueScenarioIndex]
    estimated <- sapply(fitIndex,
        # require a fitted model with method != 'none'
        # findarg new in 2.8.3 - see utility.R
        function(x) {
            ((is.logical(object$fit) && object$fit) || (object$fit %in% "multifit")) && 
                (findarg(object$fit.args, 'method', x, 'default') != 'none')
            
        }
    )
    if (missing(true)) {
        splitScenarios <- split(object$scenarios, object$scenarios$scenario)
        trueD <- sapply(splitScenarios, weighted, param='D')   ## vector length = number of scenarios
        
    if (object$outputtype == 'regionN') {
        true <- trueD
        true <- true * attr(object, 'regionsize')  ## for each scenario
    }
    else if (object$outputtype %in% c('predicted','derived')){
        if (parameter == 'D')
            true <- trueD
        else {
            if (parameter %in% names(object$scenarios))  # 2017-11-03
                true <- sapply(splitScenarios, weighted, parameter)
            else
                true <- NA
        }
    }
    else true <- NA
}
else if (length(true) != length(object$output))
    stop ("specify one 'true' value for each scenario")
    object$output <- mapply(extractfn,
                            object$output,
                            true,
                            estimated,
                            SIMPLIFY = FALSE)
    object$outputtype <- 'numeric'
    class(object) <- c('selectedstatistics', 'secrdesign', 'list')
    attr(object, 'parameter') <- parameter
    object
}
