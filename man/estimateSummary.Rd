\name{estimateSummary}
\alias{estimateArray}
\alias{estimateSummary}

\title{Direct summary of estimate tables}

\description{

An alternative approach to summarising output from \code{run.scenarios}
(cf \code{\link{summary.estimatetables}}).

}

\usage{

estimateArray(object)

estimateSummary(object, parameter = "D", statistics = c("n", "estimate", 
    "SE.estimate", "RB", "seRB", "RSE", "RMSE", "COV"), true,
    validrange = c(0, Inf), checkfields = c('estimate','SE.estimate'),
    format = c('list', 'data.frame'))

}

\arguments{
  \item{object}{secrdesign object of class "estimatetables"}
  \item{parameter}{character name of parameter (row in estimate table)}
  \item{statistics}{character choice of outputs}
  \item{true}{numeric vector of true values, one per scenario and group}
  \item{validrange}{numeric allowed for estimates or other checkfields}
  \item{checkfields}{character choice of columns in each estimate table that 
  will be checked against validrange}
  \item{format}{character choice of output}
}

\details{

When 'predict(fittedmodel)' in \code{\link{run.scenarios}} generates more than 
one estimate table (i.e. when the model uses groups, mixture classes or 
multiple sessions), the default extract function retains only the first. 
This is often OK, but it can be frustrating if you care about 
group- or session-specific estimates.

The alternative is to use 'predict' as the run.scenarios extractfn, which 
retains all estimate tables. This requires a different function for 
summarisation; \code{estimateSummary} will suffice for many purposes.

\code{estimateSummary} internally calls \code{estimateArray} to pre-process 
the output from run.scenarios.

The code should be examined for the precise definition of each statistic.

True parameter values are required for RB, RMSE and COV, and these are computed 
even if later dropped from the output. If provided, the argument \code{true} 
should have length equal to the number of parameter tables in each replicate, 
i.e. (number of scenarios) * (number of groups), ordered by scenario. 
Otherwise, true values will be taken from rows of the data frame object$scenarios.

Replicates are rejected (set to NA) if any \code{checkfields} 
falls outside \code{validrange}.

Output statistics `n', `estimate', `SE.estimate', `RB', and `RSE' are means 
across replicates.

}

\value{

\code{estimateArray} --
array with dimensions (Parameter, statistic, Group, Scenario, Replicate)\\

\code{estimateSummary} --

list of matrices (group x scenario), one for each statistic

}

\note{

These functions are experimental in version 2.8.1 and may be changed.

}

\seealso{

\code{\link{run.scenarios}},
\code{\link{summary.estimatetables}}

}

\examples{

\donttest{

scen8 <- make.scenarios (D = 8, g0 = 0.3, sigma = 30, 
    noccasions = c(4,8), groups = c('F','M'))
    
# replace density and sigma values of males to define more interesting 
# 2-scenario, 2-group simulation
male <- scen8$group == 'M'
scen8$D[male] <- 4
scen8$sigma[male] <- 40
details <- scen8[,c(1,2,4,6:8)]
details

grid <- make.grid(8, 8, spacing = 30)
mask <- make.mask(grid, buffer = 160, type = 'trapbuffer')

old <- options(digits = 3)

#--------------------------------------------------------------------------
# run a few simulations

sims <- run.scenarios(10, scen8, trapset = grid, fit = TRUE, 
    fit.args = list(model = list(D~g, g0~1, sigma~g), groups = 'group'),
    extractfn = predict, maskset = mask, ncores = 2)

# default summary uses true = c(8,4,8,4)  (cf details)
tab <- estimateSummary(sims, 'D',  c("n", "estimate", "RB", "seRB"))

# re-format list output by scenario --
data.frame(lapply(tab,t))

# format by scenario and group, and pre-pend details --
tab <- estimateSummary(sims, 'D',  c("n", "estimate", "RB", "seRB"), 
    format = 'data.frame')
data.frame(details, tab)

#--------------------------------------------------------------------------
# try with default extractfn (single table per replicate, despite groups)
sims2 <- run.scenarios(10, scen8, trapset = grid, fit = TRUE, 
     maskset = mask, ncores = 2)

# Fails with "Error in estimateSummary(sims2, "D") : incongruent 'true'""
# tab2 <- estimateSummary(sims2, 'D')
# OK if manually provide scenario-specific true density
tab2 <- estimateSummary(sims2, 'D', true = c(12,12))

# reformat by scenario
estimateSummary(sims2, 'D', true = c(12,12), format = 'data.frame')

# compare
summary(sims2)
#--------------------------------------------------------------------------

# multiple estimate tables also arise from multi-session simulations
# argument 'true' must be specified manually
# interpret with care: sessions are (probably) not independent
# this example uses the previous grid and mask
scen9 <- make.scenarios (D = 8, g0 = 0.3, sigma = 30, noccasions = 5)
poparg <- list(nsessions = 3, details = list(lambda = 1.2))  # for sim.popn
detarg <- list(renumber = FALSE)                             # for sim.capthist
fitarg <- list(model = D~Session)                            # for secr.fit
sims3 <- run.scenarios(5, scen9, trapset = grid, fit = TRUE, 
    maskset = mask, ncores = 2, pop.args = poparg, det.args = detarg,
    fit.args = fitarg, extractfn = predict)
estimateSummary(sims3, parameter = 'D', format = 'data.frame', true = 8 * 1.2^(0:2))

# extractfn = coef usually results in a single estimate table,
# so the usual summary method is usually sufficent.
# For completeness we show that estimateSummary can also be used.
# Coefficients are often negative, so relative values may be meaningless
sims4 <- run.scenarios(5, scen9, trapset = grid, fit = TRUE, 
    maskset = mask, ncores = 2, pop.args = poparg, det.args = detarg,
    fit.args = fitarg, extractfn = coef)
estimateSummary(sims4, parameter = 'D', c("n", "estimate", "SE.estimate", "RMSE", "COV"), 
    format = 'data.frame', true = log(8))

options(old)

}

}
