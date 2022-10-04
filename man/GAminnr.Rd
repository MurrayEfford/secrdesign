\name{GAminnr}
\alias{GAminnr}

\title{

SECR detector placement by min(n,r) optimization.

}

\description{

Implements the approach of Durbach et al. (2021) for optimization of detector 
placement using a genetic algorithm. 

}

\usage{

GAminnr(mask, alltraps, ntraps, detectpar, noccasions = 1, 
    detectfn = c("HHN", "HHR", "HEX", "HAN", "HCG"), D = 1, 
    criterion = 4, pen_gridsigma = 2, pen_wt = 0, 
    seed = NULL, ...)
    
}

\arguments{

  \item{mask}{mask object}
  \item{alltraps}{traps object with all possible trap locations}
  \item{ntraps}{number of required trap locations }
  \item{detectpar}{list values of detection parameters lambd0, sigma etc.}
  \item{detectfn}{integer code or character string for shape of detection function - see \link{detectfn}}
  \item{noccasions}{integer number of sampling occasions}
  \item{D}{numeric density animals/ha}
  \item{criterion}{integer criterion to use (1 = En, 2 = Er, 3 = Em, 4 = min(En,Er))}
  \item{pen_gridsigma}{numeric penalty for poor trap spacing }
  \item{pen_wt}{numeric penalty for spacings that are too different from regular grid}
  \item{seed}{set a random seed for reproducibility of GA}
  \item{\dots}{other arguments passed to \code{\link{kofnGA}}}

}

\details{

\code{detectpar} is a named list with values of the detection parameters for the chosen detectfn. Usually this will be just lambda0 (baseline hazard of detection) and sigma (spatial scale of detection).

If pen_wt>0 then penalties for poor trap spacing are based on comparing the min(n,r) grid to what you would get from a regular grid with spacing sigma * pen_gridsigma. [unclear]

The genetic algorithm is provided by function \code{kofnGA} from package \pkg{kofnGA} (Wolters 2015). The first three arguments of \code{kofnGA} (i.e., n, k, OF) are set by \code{GAminnr}. Others may be adjusted by the user via the \dots argument. Specifically,

\tabular{lll}{

Argument \tab Default \tab Description \cr
ngen \tab 500 \tab number of generations to run \cr

popsize \tab 200 \tab size of the population; equivalently, the number of offspring produced each generation \cr

mutprob \tab 0.01 \tab mutation rate \cr

verbose \tab 0 \tab integer controlling the display of progress during search. If a positive value, then the iteration number and best objective function value are displayed at the console every 'verbose' generations. Otherwise nothing is displayed. The default gives no display. \cr

cluster \tab NULL \tab number of parallel cores or a prebuilt parallel cluster \cr

}

The default for ngen may (or may not) be too large for routine use. Durbach et al. (2021) used ngen = 50, popsize = 1000 and 
mutprob = 0.01. 

}

\value{

GAminnr -- A list with components

\item{mask}{saved input}
\item{alltraps}{saved input}
\item{detectpar}{saved input}
\item{noccasions}{saved input}
\item{detectfn}{saved input}
\item{D}{saved input}
\item{des}{kofnGA() output object}
\item{optimaltraps}{traps object with optimized layout}
\item{optimalenrms}{E(n), E(r), E(m) evaluated with optimized layout}

}

\references{

Durbach, I., Borchers, D., Sutherland, C. and Sharma, K. (2021) Fast, flexible alternatives to regular grid designs for spatial capture--recapture. \emph{Methods in Ecology and Evolution} \bold{12}, 298--310. DOI 10.1111/2041-210X.13517

Wolters, M. A. (2015) A genetic algorithm for selection of fixed-size
subsets with application to design problems. \emph{Journal of Statistical
Software, Code Snippets}, \bold{68}, 1--18. DOI 10.18637/jss.v068.c01

}

\author{

Ian Durbach. Adapted for \pkg{secrdesign} by Murray Efford.

}

\seealso{

\code{\link{Enrm}}

}

\examples{

\donttest{

}

}
