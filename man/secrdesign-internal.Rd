\name{Internal}

\alias{penfn}

\title{Internal Functions}

\description{ 

Functions that may be called internally by \pkg{secrdesign}. 
These are exported and may be called separately for testing.

}

\usage{

penfn(traps, sigma)

}

\arguments{

  \item{traps}{secr trapsobject}
  \item{sigma}{numeric sparial scale parameter}
  
}

\details{

  \code{penfn} is the default pen_fn used by \code{\link{GAminnr}} When called with a non-null penalty argument.

}

\value{

penfn -- a numeric vector with the number of trap pairs separated by 2.5-3.5 sigma and 3.5-4.5 sigma.

}

\seealso{

\code{\link{GAminnr}}, 

}

\references{

Durbach, I., Borchers, D., Sutherland, C. and Sharma, K. (2021) Fast, flexible alternatives to regular grid designs for spatial capture--recapture. \emph{Methods in Ecology and Evolution} \bold{12}, 298--310. DOI 10.1111/2041-210X.13517

}

\examples{

penfn(traps(captdata), 20)

}