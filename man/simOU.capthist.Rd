\name{simOU.capthist}
\alias{simOU}
\alias{simOU.capthist}

\title{Simulate Detections For Ornstein-Uhlenbeck Movement}

\description{

The usual SECR model ignores the sequential locations of an individual within 
its home range. Movement models predict serial correlation of detections in space. 
The Ornstein-Uhlenbeck (OU) model is a convenient example that over long 
durations leads to a bivariate normal distribution of locations.

Movements of a single animal according to the OU model are simulated in discrete time with \code{simOU}. 

Detections of a population of individuals with pre-defined activity centres are 
simulated with \code{simOU.capthist}. Detection happens when the location of an 
individual at time \eqn{t} (occasion \eqn{t}) is within distance epsilon \eqn{\epsilon} of a detector.

}

\usage{

simOU(xy, tau, sigma, noccasions, start = NULL)

simOU.capthist(traps, popn, detectpar, noccasions, epsilon, 
    savepopn = FALSE, savepath = FALSE, ...)

}

\arguments{

  \item{xy}{numeric vector of x,y coordinates for one animal}
  \item{tau}{numeric parameter of serial correlation = \eqn{1 \ \beta}}
  \item{sigma}{numeric spatial scale parameter}
  \item{noccasions}{integer number of time steps}
  \item{start}{numeric vector of x,y coordinates for initial location (optional)}

  \item{traps}{secr traps object}
  \item{popn}{secr popn object or a 2-column matrix of x-y coordinates of activity centres}
  \item{detectpar}{list with values of detection parameters tau and sigma}
  \item{epsilon}{numeric radius within which individual detected with certainty}
  \item{savepopn}{logical; if TRUE the population is saved as an attribute}
  \item{savepath}{logical; if TRUE the movement paths are saved as an attribute}
  \item{\dots}{other arguments passed to \code{\link[secr]{reduce.capthist}}}
}

\details{

The first location for \code{simOU} by default is drawn at random from the asymptotic distribution.

If 'detectpar' does not include a component 'tau' then the value of tau will be 
taken from detectpar$lambda0.

In a final step, `simOU.capthist' uses the reduce method for capthist objects to 
coerce the simulated capthist object to the detector type of the traps argument. 

The \dots argument may be used to pass 'by = "ALL"' to reduce.capthist, 
collapsing the binary data for a single detector on noccasions to a count. 
This requires that the input traps object has detector type 'count'.

}

\value{

simOU - matrix of locations dim = c(noccasions, 2)

simOU.capthist - single-session capthist object for \pkg{secr}

}

\seealso{

\code{\link[secr]{sim.capthist}}

}

\examples{

# one animal
locs <- simOU(c(0,0), 20, 1, 100)
par(pty = 's')
plot(locs, type = 'o', xlim = c(-2.5,2.5), ylim = c(-2.5,2.5))
points(0,0, pch = 16, col = 'red')

\donttest{

# simulate some capture data
set.seed(123)
grid <- make.grid(8, 8, spacing = 2)
pop <- sim.popn(D = 1000, core = grid, buffer = 4)
ch <- simOU.capthist(grid, pop, detectpar=list(tau = 50, sigma = 1), 
   noccasions = 100, epsilon = 0.25, savepath = TRUE)

# plot simulated capthist with overlay of movements and AC
plot(ch, rad = 0.1, tracks = TRUE, varycol = FALSE, border = 4)
sapply(attr(ch, 'path'), lines, col = 'grey')
plot(pop, add = TRUE, pch = 16, cex = 0.6)

# fit a model
fit <- secr.fit(capthist = ch, buffer = 8, detectfn = 14, trace = FALSE)
predict(fit)

}

}