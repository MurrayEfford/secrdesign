\name{defaultextractfn}
\alias{defaultextractfn}

\title{Default Extract Function for run.scenarios}

\description{

The function applied by \code{\link{run.scenarios}} to the result of each
replicate when \code{extractfn} is not provided. It may be called from a
user-written extract function that modifies or adds to the default output.

}

\usage{

defaultextractfn(x, ...)

}

\arguments{

  \item{x}{object to summarise: a fitted model (class `secr' or `ipsecr'),
  a list of fitted models (class `secrlist'), a capthist object, or a
  `try-error' from a failed fit}
  \item{...}{other arguments, including \code{keepdensity} (see Details)}

}

\details{

The result depends on the class of \code{x}.

\tabular{ll}{
fitted model \tab data.frame from \code{predict(x)} with attributes \code{counts}, \code{start}, \code{method}, \code{proctime}, \code{chat}, \code{code} and \code{convergence}\cr
secrlist \tab named list of the above, one per model (see below)\cr
capthist \tab vector of counts of the raw data \cr
try-error \tab data.frame with no rows or columns \cr
}

If \code{keepdensity = TRUE} is passed to \code{run.scenarios} then the
simulated density surface (attribute \code{D} of the capthist object, if present) is
saved as attribute \code{D}.

With \code{fit = "multifit"}, \code{run.scenarios} fits several models to each
simulated dataset and \code{x} is a \code{secrlist}. Models are named by the
names of the nested \code{fit.args} list if provided; otherwise by their
formulae if these differ, or \code{model1}, \code{model2}, etc.
A user-written \code{extractfn} for \code{fit = "multifit"}
should return a named list of data.frames of the same form
(one per model) if the results are to be summarised with
\code{\link{estimateSummary}} or \code{\link{countSummary}}.
If a fit fails in a replicate, \code{x} is a plain list that includes a `try-error'
and the corresponding component of the result is an empty data.frame.

}

\value{

See Details.

}

\seealso{

\code{\link{run.scenarios}}, \code{\link{estimateSummary}}

}

\examples{

\dontrun{

# two models fitted to each simulated dataset
fitargs <- list(
    null = list(model = list(D~1, g0~1, sigma~1)),
    bhv  = list(model = list(D~1, g0~b, sigma~1)))
# fit.args is a nested list: one list of models for each fitindex
sims <- run.scenarios(100, scen, trps, msk, fit = "multifit", fit.args = list(fitargs))
estimateSummary(sims, "D")                 # all models, with a 'model' column
estimateSummary(sims, "D", model = "bhv")

}

}
