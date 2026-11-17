\name{Etcc}
\alias{Etcc}
\alias{RSEtcc}

\title{Expected Counts and Precision for Telemetry-Calibrated Cue Counts}

\description{

Compute the expected counts of cues, and an approximate RSE of the density estimate,
for a study design that combines cue counts at `count' detectors (cameras) with a
sample of radio-collared animals (telemetry type `marking' in \pkg{secr};
see \code{\link[secr]{addTelemetry}}). \code{Etcc} is the counterpart of
\code{\link{Enrm}} for this method: it takes the parameters of detection, the density,
the detector array, the number of occasions and the number of collars, and needs no
simulation. \code{RSEtcc} has the same inputs and returns the expected counts with
approximations to the relative standard error of the density estimate, as
\code{\link{minnrRSE}} does for \code{Enrm}.

}

\usage{

Etcc(D, lambda0, sigma, traps, mask, noccasions = 5, ncollar = 20,
    deployment = c("pdot", "uniform"), q = 1, nfix = 100,
    Ndist = c("poisson", "fixed"), warnseen = 3, detectfn = c("HHN", "HEX"))

RSEtcc(D, lambda0, sigma, traps, mask, noccasions = 5, ncollar = 20,
    deployment = c("pdot", "uniform"), q = 1, nfix = 100,
    Ndist = c("poisson", "fixed"), calib = c("mc", "fisher", "poisson"),
    kcal = 1, warnseen = 3, nsimcollar = 2000, detectfn = c("HHN", "HEX"))

}

\arguments{
  \item{D}{population density animals / hectare; scalar, or a vector of length \code{nrow(mask)} for a single setting}
  \item{lambda0}{numeric intercept of the hazard detection function}
  \item{sigma}{numeric spatial scale of detection (same length units as \code{traps})}
  \item{traps}{\code{\link[secr]{traps}} object with detector type `count'}
  \item{mask}{\code{\link[secr]{mask}} object (the buffered habitat region)}
  \item{noccasions}{integer number of sighting occasions}
  \item{ncollar}{integer number of collared animals}
  \item{deployment}{character; how collars are distributed (see Details)}
  \item{q}{numeric probability that a sighting of a collared animal is identified}
  \item{nfix}{integer number of telemetry fixes per collared animal}
  \item{Ndist}{character; distribution of the number of animals (see Details)}
  \item{warnseen}{numeric; \code{warn} is TRUE if fewer collared animals than this are expected to be detected}
  \item{calib}{character; how the contribution of the collared sample is computed (see Details)}
  \item{kcal}{numeric multiplier of 1/R used when \code{calib = "poisson"}}
  \item{nsimcollar}{integer number of random sets of collars when \code{calib = "mc"}}
  \item{detectfn}{character; shape of the hazard of detection: half-normal \code{"HHN"} or exponential \code{"HEX"} (see Details)}
}

\details{

The arguments \code{D}, \code{lambda0}, \code{sigma}, \code{noccasions}, \code{ncollar}, \code{q} and \code{nfix} are recycled to a common length; each element defines one setting and gives one row of the result. Use \code{\link{expand.grid}} to combine settings. Only the hazard half-normal (\code{"HHN"}) and hazard exponential (\code{"HEX"}) detection functions are implemented (as for telemetrytype 'marking' in \pkg{secr}), and detectors are assumed to be used on every occasion.

\subsection{Expected counts}{

Let \eqn{L(x) = S \sum_k \lambda_0 \exp(-d_k(x)^2/2\sigma^2)} (\code{"HHN"}; \eqn{\exp(-d_k(x)/\sigma)} for \code{"HEX"}) be the expected number of cues in \eqn{S} occasions from an animal with activity centre at \eqn{x}. With \eqn{a} the cell area of the mask, define \eqn{J_1 = a \sum_x D(x) L(x)} (all cues), \eqn{J_2 = a \sum_x D(x) L(x)^2} and \eqn{N = a \sum_x D(x)}. The expected number of animals detected at least once is \eqn{E(n) = a \sum_x D(x) (1 - \exp(-L(x)))}, as for any SECR design (see \code{\link{Enrm}}).

The collared animals' activity centres are distributed with probability \eqn{w(x)} over the mask: proportional to \eqn{D(x)} if \code{deployment = "uniform"} (widespread), or to \eqn{D(x)} times the probability of detection in one occasion if \code{deployment = "pdot"} (centred; an analogue of selecting animals that are detected). The expected number of cues from the collared animals is \eqn{C = m \sum_x w(x) L(x)} for \eqn{m} collars. Then

\itemize{
\item \code{U} = \eqn{J_1 - C} expected cues of unmarked animals (the total \eqn{T^u}),
\item \code{R} = \eqn{qC} expected identified cues of collared animals, and \code{Tm} = \eqn{(1-q)C} expected unidentified cues,
\item \code{seen} expected number of collared animals detected at least once,
\item \code{collarshare} = \eqn{C/J_1}, the share of all cues that come from collared animals.
}

}

\subsection{Overdispersion}{

The total unmarked count is a compound Poisson variable; its variance to mean ratio is
\deqn{\hat c = 1 + (J_2 + V_s)/\mu}{c = 1 + (J2 + Vs) / mu}
where \eqn{\mu} = \code{U}. \eqn{V_s = m \sum_x w(x) |\nabla L|^2 v/n_{fix}} approximates the variance in the collared animals' cues caused by uncertainty in their activity centres; \eqn{v} is the per-axis variance of an activity centre estimated from one fix (\eqn{\sigma^2} for \code{"HHN"}, \eqn{2\sigma^2} for \code{"HEX"}). With \code{Ndist = "fixed"}, \eqn{J_2 - J_1^2/N} replaces \eqn{J_2}. This is the expression used by \code{secr.fit} with \code{details = list(chatmethod = "analytic")} (\code{chat}).

}

\subsection{Predicted RSE}{

\code{RSEtcc} gives two approximations. The coarse approximation, in the spirit of \code{minnrRSE}, adds the relative variances from the animals detected and from the collar cues:
\deqn{\mathrm{RSE} \approx \sqrt{1/E(n) + 1/C}.}{RSE = sqrt(1/E(n) + 1/C).}
This is \code{RSEcoarse}. The closer approximation is based on density being estimated as total cues divided by the expected cues per animal, so two components of variance are added:
\deqn{\mathrm{RSE}^2 = \hat c \mu/J_1^2 + \mathrm{CV}^2(\hat\Lambda).}{RSE^2 = c mu / J1^2 + CV^2(Lambda-hat).}
The first term (\code{RSEunmarked}) comes from the count of unmarked cues; it changes little with the number of collars, because removing collared animals lowers the variance of the unmarked count and the mean together. The second (\code{RSEcalib}) comes from estimating the cue rate of an animal from the collared sample, and decreases roughly as 1/C. Overall \code{RSE} is the square root of the sum.

The method for the second term is set by \code{calib}.

\describe{
\item{"mc"}{(default) The information in the collared sample depends on where the collars happen to be. Sets of \eqn{m} activity centres are drawn from \eqn{w(x)}, the number of collar cues is drawn from a Poisson distribution given the centres, and \eqn{\mathrm{CV}^2} is the mean of 1/cues over draws with at least one cue. This allows for the sets of collars that are far from the cameras and carry little information, which is important when deployment is widespread or collars are few. The proportion of draws with no cues is returned as \code{pnocues}; there the model could not be fitted. Results vary slightly among calls unless \code{set.seed} is used.}
\item{"fisher"}{The delta-method variance of the mean expected cues per animal from the expected Fisher information for \eqn{\lambda_0} and \eqn{\sigma} (collared animals' counts, with activity centres treated as known, plus the telemetry fixes). For the designs examined this equals 1/(collar cues) closely, because the telemetry determines \eqn{\sigma}; it does not depend on \code{q}.}
\item{"poisson"}{\code{kcal} / \code{R}, using the expected value of R rather than its distribution.}
}

}

}

\value{

A data frame with one row for each setting. \code{Etcc} returns the columns

  \item{D, lambda0, sigma, noccasions, ncollar, deployment, q}{settings}
  \item{N}{expected number of animals in the mask}
  \item{En}{expected number of animals detected at least once}
  \item{J1}{expected total number of cues from all animals}
  \item{U}{expected number of cues of unmarked animals}
  \item{R}{expected number of identified cues of collared animals}
  \item{Tm}{expected number of unidentified cues of collared animals}
  \item{seen}{expected number of collared animals detected at least once}
  \item{collarshare}{share of cues from collared animals}
  \item{chat}{variance to mean ratio of the unmarked total}
  \item{warn}{TRUE if \code{seen < warnseen}}

\code{RSEtcc} returns these columns and in addition

  \item{RSEcoarse}{coarse approximation to the RSE of the density estimate}
  \item{RSEunmarked}{component of RSE from the unmarked count}
  \item{RSEcalib}{component of RSE from the collared sample}
  \item{RSE}{predicted RSE of the density estimate}
  \item{pnocues}{proportion of random sets of collars with no cues (\code{calib = "mc"} only)}

}

\section{Limitations}{

The predictions have been compared with simulated density estimates (secr 5.5.2, \code{telemetrytype = "marking"}) only for homogeneous density, one class of animals, \eqn{q = 1}, all collars active for the whole survey, and a 9 x 9 grid of cameras; \code{RSE} was mostly within 10--20\% of the simulated RSE and tended to be high for centred deployment, and \code{RSEcoarse} was between about 0.7 and 1.1 of it, lowest for few widespread collars. Neither allows for clumped activity centres or other spatial variation in density beyond what is supplied in \code{D}, for overdispersion in the cue counts other than that from animals, for heterogeneity among animals, or for the bias that results from few collars (the `ratio effect'). Predictions for fewer than about 5 expected collar detections are optimistic (see \code{warn}); estimates may also fail.

The RSE is for the density of the population over \code{mask}, assuming collared animals are representative of the population in detectability.

}

\references{

Efford, M. G., and Boulanger, J. (2019) Fast evaluation of study designs for
spatially explicit capture--recapture.
\emph{Methods in Ecology and Evolution}, \bold{10}, 1529--1535. \doi{10.1111/2041-210X.13239}

}

\seealso{

\code{\link{Enrm}}, \code{\link{minnrRSE}}, \code{\link{Echat}}, \code{\link{Lambda}}, \code{\link[secr]{addTelemetry}}, \code{\link[secr]{details}}

}

\examples{

tr <- make.grid(nx = 9, ny = 9, spacing = 11.25, detector = "count", origin = c(60,60))
msk <- make.mask(tr, buffer = 60, type = "traprect", spacing = 5)

# expected counts
Etcc(D = 53, lambda0 = 0.54, sigma = 8, traps = tr, mask = msk, ncollar = 20)

# predicted RSE vs number of collars, widespread deployment
set.seed(1)
RSEtcc(D = 53, lambda0 = 0.54, sigma = 8, traps = tr, mask = msk,
    ncollar = c(5, 10, 20, 40), deployment = "uniform")[, c("ncollar", "En", "seen", 
    "RSEcoarse", "RSE", "warn")]

\donttest{
n <- 3:60
rse <- sapply(c("pdot", "uniform"), function(dep)
    RSEtcc(D = 53, lambda0 = 0.54, sigma = 8, traps = tr, mask = msk, ncollar = n,
         deployment = dep)$RSE)
matplot(n, rse, type = "l", lty = 1, xlab = "Number of collars", ylab = "RSE(D-hat)")
legend("topright", c("centred", "widespread"), col = 1:2, lty = 1)
}

}
