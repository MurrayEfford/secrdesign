\name{Echat}
\alias{Echat}

\title{Expected Variance Inflation Factor for Sighting Data}

\description{

Compute the expected variance inflation factor (c-hat) of the pooled total of 
sightings of unmarked animals, from the parameters of detection, the density, 
the detector array and the number of occasions, without simulation. This is the 
design-stage counterpart of the analytic c-hat that \pkg{secr} computes at the 
fitted parameters (\code{details = list(chatmethod = "analytic")} in 
\code{\link[secr]{secr.fit}}), and a companion to \code{\link{Enrm}}.

}

\usage{

Echat(D, traps, mask, detectpar, noccasions,
    detectfn = c("HHN", "HHR", "HEX", "HN", "HR", "EX"),
    Ndist = c("poisson", "fixed"))

}

\arguments{
  \item{D}{population density animals / hectare; scalar, or a vector of length \code{nrow(mask)}}
  \item{traps}{\code{\link[secr]{traps}} object with detector type `count' or `proximity'}
  \item{mask}{\code{\link[secr]{mask}} object}
  \item{detectpar}{list of detection parameters (as for \code{\link{Enrm}})}
  \item{noccasions}{integer number of sighting occasions}
  \item{detectfn}{character or integer; detection function}
  \item{Ndist}{character; distribution of the number of animals}
}

\details{

The model is that of the `sighting-only' models of \pkg{secr}: all animals are unmarked, activity centres are independent with distribution proportional to \eqn{D(x)}, and the number of animals is Poisson or fixed. All detectors are assumed to be used on every occasion (as in \code{Enrm}).

Cues of one animal fall at several detectors and occasions, so the total \eqn{T} is not a Poisson variable. The variance inflation factor is \eqn{\hat c = \mathrm{Var}(T) / \mathrm{E}(T)} for count detectors, and for proximity detectors the ratio of \eqn{\mathrm{Var}(T)} to the variance if the cells were independent Bernoulli variables with a common probability, where \eqn{T} is the number of detector-occasions with a detection (`positive cells').

\emph{Count detectors.} Let \eqn{L(x)} be the expected number of cues from an animal centred at \eqn{x} (summed over detectors and occasions), \eqn{a} the cell area of the mask, \eqn{J_1 = a \sum_x D(x) L(x)} and \eqn{J_2 = a \sum_x D(x) L(x)^2}. Then \eqn{\mathrm{E}(T) = J_1} and \eqn{\mathrm{Var}(T) = J_1 + J_2} for Poisson \eqn{N}, and \eqn{J_1 + J_2 - J_1^2/N} for fixed \eqn{N}.

\emph{Proximity detectors.} A cell \eqn{(k,s)} is positive with probability \eqn{1 - \exp(-X)}, where \eqn{X} is the summed hazard of the animals present. With \eqn{G} the probability generating function of the number of animals, \eqn{u(x) \propto D(x)}, \eqn{h_k(x)} the hazard at detector \eqn{k}, and
\deqn{M_{1k} = G\{\sum_x u(x) \exp(-h_k(x))\}, \quad M_{2kk'} = G\{\sum_x u(x) \exp(-h_k(x) - h_{k'}(x))\},}{M1[k] = G(sum_x u(x) exp(-h[k](x))), M2[k,k'] = G(sum_x u(x) exp(-h[k](x) - h[k'](x))),}
the mean is \eqn{S \sum_k (1 - M_{1k})} and the variance is \eqn{S \sum_k (M_{1k} - M_{2kk}) + S^2 \sum_{kk'} M_{2kk'} - S^2 (\sum_k M_{1k})^2}, with \eqn{S} the number of occasions.

The result agrees with the value computed by \pkg{secr} for fitted models, and with simulation.

}

\value{

Named numeric vector with elements

  \item{ET}{expected value of the total: the number of cues (count) or of positive cells (proximity)}
  \item{VT}{variance of the total}
  \item{chat}{variance inflation factor}

}

\note{

\code{Echat} does not allow for marked animals (cf. \code{\link{Etcc}}, which removes the collared animals' expected cues for telemetry-calibrated cue counts), for varying usage, or for `multi' detectors.

}

\references{

Efford, M. G. and Hunter, C. M. (2018) Spatial capture--mark--resight estimation of animal population density. \emph{Biometrics} \bold{74}, 411--420.

}

\seealso{

\code{\link{Enrm}}, \code{\link{Etcc}}, \code{\link[secr]{details}}

}

\examples{

tr <- make.grid(nx = 6, ny = 6, spacing = 15, detector = "count")
msk <- make.mask(tr, buffer = 45, type = "traprect", spacing = 7.5)
dp <- list(lambda0 = 0.4, sigma = 15)

# count detectors
Echat(D = 20, traps = tr, mask = msk, detectpar = dp, noccasions = 5)

# binary proximity detectors, same parameters
detector(tr) <- "proximity"
Echat(D = 20, traps = tr, mask = msk, detectpar = dp, noccasions = 5)

# fixed number of animals
Echat(D = 20, traps = tr, mask = msk, detectpar = dp, noccasions = 5, Ndist = "fixed")

}
