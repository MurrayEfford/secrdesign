\name{rbind.estimatetables}
\alias{rbind.estimatetables}
\alias{c.estimatetables}

\title{

Combine Simulation Output

}

\description{

Methods to combine estimate tables from separate executions of \code{\link{run.scenarios}}.

}

\usage{

\method{rbind}{estimatetables}(..., deparse.level = 1)

\method{c}{estimatetables}(...)

}

\arguments{
  \item{\dots}{estimatetables output from \code{\link{run.scenarios}}}
  \item{deparse.level}{not used (required by genberic method \code{\link{rbind}})}
}

\details{
\code{rbind} assumes all inputs used exactly the same scenarios. Replicate estimate tables are combined across executions for each scenario in turn. This is useful to increase the number of replicates. The `scenarios' component remains unchanged.

\code{c} combines outputs from \code{run.scenarios} that may differ in their scenarios. The `output' component of the result is a concatenation of the output lists in the input. The `scenarios' component of the result comprises the input scenarios stacked with \code{\link{rbind.data.frame}}.

The compatibility of the inputs is checked, but the checks are not exhaustive, and users should be wary.

}

\value{

`estimatetables' object combining the inputs

}

\seealso{

\code{\link{make.scenarios}}
\code{\link{run.scenarios}}

}

\examples{
}
