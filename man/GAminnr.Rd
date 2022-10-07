\name{GAminnr}
\alias{GAminnr}

\title{

SECR detector placement by min(n,r) optimization.

}

\description{

Implements the approach of Durbach et al. (2021) for optimization of detector 
placement using a genetic algorithm to maximize the lesser 
of E(\eqn{n}) and E(\eqn{r}), where \eqn{n} is the number of distinct individuals 
and \eqn{r} is the total number of recaptures. This criterion predicts the 
relative standard error of the density estimate (Efford and Boulanger 2019).

}

\usage{

GAminnr(mask, alltraps, ntraps, detectpar, noccasions, 
    detectfn = c("HHN", "HHR", "HEX", "HAN", "HCG"), D = NULL, 
    penalty = NULL, seed = NULL, ...)
    
}

\arguments{

  \item{mask}{mask object}
  \item{alltraps}{traps object with all possible trap locations}
  \item{ntraps}{number of required trap locations }
  \item{detectpar}{list values of detection parameters lambd0, sigma etc.}
  \item{detectfn}{integer code or character string for shape of detection function - see \link{detectfn}}
  \item{noccasions}{integer number of sampling occasions}
  \item{D}{numeric density animals per hectare (0.01 km^2)}
  \item{penalty}{list defining penalty for layout in relation to reference grid (optional)}
  \item{seed}{set a random seed for reproducibility of GA}
  \item{\dots}{other arguments passed to \code{\link{kofnGA}}}

}

\details{

\code{detectpar} is a named list with values of the detection parameters for the chosen detectfn. Usually this will be just lambda0 (baseline hazard of detection) and sigma (spatial scale of detection).

The genetic algorithm is provided by function \code{kofnGA} from package \pkg{kofnGA} (Wolters 2015). The first three arguments of \code{kofnGA} (i.e., n, k, OF) are set by \code{GAminnr}. Others may be adjusted by the user via the \dots argument. Specifically,

\tabular{lll}{

Argument \tab Default \tab Description \cr
ngen \tab 500 \tab number of generations to run \cr

popsize \tab 200 \tab size of the population; equivalently, the number of offspring produced each generation \cr

mutprob \tab 0.01 \tab mutation rate \cr

verbose \tab 0 \tab integer controlling the display of progress during search. If a positive value, then the iteration number and best objective function value are displayed at the console every 'verbose' generations. Otherwise nothing is displayed. The default gives no display. \cr

cluster \tab NULL \tab number of parallel cores or a prebuilt parallel cluster \cr

}

The default for ngen may (or may not) be larger than is needed for routine use. Durbach et al. (2021) used ngen = 50, popsize = 1000 and mutprob = 0.01. 

Density D may be a scalar or a vector of length equal to the number of mask cells. No value need be specified if the sole aim is to optimize trap placement, but D is required for predictions of E(\eqn{n}) and E(\eqn{r}).

Pathological trap layouts (sensu Efford and Boulanger 2019) may be avoided by adding a penalty to the objective.
No penalty is applied by default. To apply a penalty, \code{penalty} should be a list with named components pen_wt>0 and pen_gridsigma). If a penalty is applied, the default compares the number of trap pairs with close spacing (2.5-3.5 sigma, 3.5-4.5 sigma) to the number in a compact sample from a regular grid with spacing sigma * pen_gridsigma (see internal functions \code{\link{GApenfn}} and \code{\link{compactSample}} and the vignette).  An alternative penalty function may be supplied as component `pen_fn' of \code{penalty}.

}

\value{

An object of class "GAminnr" that is a list with components

\item{mask}{saved input}
\item{alltraps}{saved input}
\item{detectpar}{saved input}
\item{noccasions}{saved input}
\item{detectfn}{saved input}
\item{D}{saved input}
\item{penalty}{saved input}
\item{des}{kofnGA() output object}
\item{optimaltraps}{traps object with optimized layout}
\item{optimalenrms}{E(n), E(r), E(m) evaluated with optimized layout}

}

\references{

Durbach, I., Borchers, D., Sutherland, C. and Sharma, K. (2021) Fast, flexible 
alternatives to regular grid designs for spatial capture--recapture. 
\emph{Methods in Ecology and Evolution} \bold{12}, 298--310. DOI 10.1111/2041-210X.13517

Efford, M. G., and Boulanger, J. (2019) Fast evaluation of study designs for 
spatially explicit capture--recapture. 
\emph{Methods in Ecology and Evolution}, \bold{10}, 1529--1535. DOI: 10.1111/2041-210X.13239

Wolters, M. A. (2015) A genetic algorithm for selection of fixed-size
subsets with application to design problems. \emph{Journal of Statistical
Software, Code Snippets}, \bold{68}, 1--18. DOI 10.18637/jss.v068.c01

}

\author{

Ian Durbach and Murray Efford.

}

\seealso{

\code{\link{Enrm}}, 
\code{\link{minnrRSE}}, 
\code{\link{GApenfn}}, 
\code{\link{compactSample}}

}

\examples{

\donttest{

# an artificial example
msk <- make.mask(type = 'rectangular', spacing = 10, nx = 30, ny = 20, buffer = 0)
alltrps <- make.grid(nx = 29, ny = 19, origin = c(10,10), spacing = 10)
set.seed(123)

# 50 generations for demonstration, use more in practice
opt <- GAminnr(msk, alltrps, ntraps = 20, detectpar = list(lambda0 = 0.5, sigma = 20), 
    detectfn = 'HHN', D = 10, noccasions = 5, ngen = 50, verbose = 1)

plot(msk)
plot(opt$optimaltraps, add = TRUE)
minnrRSE(opt, distribution = 'binomial')

}

}
