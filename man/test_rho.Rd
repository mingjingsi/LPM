\name{test_rho}
\alias{test_rho}
\title{
  Relationship test among traits
}
\description{
Relationship test among traits.
}
\usage{
test_rho(bLPMfit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bLPMfit}{
     bLPM fit.
}
}
\details{
  \code{test_rho} implements the relationship test among traits.
}
\value{
Returns \code{p_value} p-value of the relationship test.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming
}
\seealso{
 \code{\link{bLPM}}, \code{\link{test_rho}}
}
\examples{
p_value_test_rho <- test_rho(bLPMfit)
}