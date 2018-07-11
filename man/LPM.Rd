\name{LPM}
\alias{LPM}
\title{
  Get the parameter estimation in LPM
}
\description{
Get the parameter estimation in LPM.
}
\usage{
LPM(bLPMfit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bLPMfit}{
     bLPM fit.
}
}
\details{
  \code{LPM} obtain the parameter estimation for LPM by using the bLPM fit.
}
\value{
List of parameter estimation.
  \item{alpha}{
     Estimation of parameter \code{alpha}.
}
  \item{beta}{
     Estimation of parameter \code{beta}.
}
  \item{R}{
     Estimation of the correlation matrix \code{R}.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming
}
\seealso{
 \code{\link{bLPM}}, \code{\link{post}}, \code{\link{test_beta}}
}
\examples{
LPMfit <- LPM(bLPMfit)
}