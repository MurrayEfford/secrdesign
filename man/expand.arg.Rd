\name{expand.arg}
\alias{expand.arg}

\title{Expand Argument List}

\description{

Generate a list of lists from vectors of argument values.

}

\usage{

expand.arg(...)

}

\arguments{

  \item{\dots}{named vectors of argument values}

}

\details{

The full 'pop.args', 'det.args' and 'fit.args' arguments of \code{\link{run.scenarios}} are lists of lists corresponding to the popindex, detindex and fitindex columns in scenarios.

\code{expand.arg} constructs such lists from all possible combinations of specified arguments, with invariant arguments appended. 

\code{expand.arg} does not work with arguments that are themselves lists, such as detectpar.

}

\value{

A list of lists.

}

\seealso{

\code{\link{run.scenarios}}

}

\examples{

## generate some arguments for sim.capthist
## noccasions is constant

expand.arg(detectfn = c('HN','EX'), noccasions = 5)

}

