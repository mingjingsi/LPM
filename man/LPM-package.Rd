\name{LPM-package}
\alias{LPM-package}
\docType{package}
\title{
LPM (Latent Probit Model)
}
\description{
This package provides functions for fitting LPM, a latent probit model to characterize relationship among complex traits using summary statistics from multiple GWASs and functional annotations.
}
\details{
\tabular{ll}{
Package: \tab LPM\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2018-07-10\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This package contains three functions \code{LPM}, \code{bLPM} and \code{bLPM_add} to fit LPM and four functions \code{post}, \code{assoc}, \code{test_rho} and \code{test_beta} to make statistical inference for risk SNPs, relationship test and hypothesis testing of annotation enrichment.
}

\author{
  \packageAuthor{LPM}
  
  Maintainer: \packageMaintainer{LPM}
}
\references{

}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
  \code{\link{LPM}}, \code{\link{bLPM}}, \code{\link{bLPM_add}}, \code{\link{post}}, \code{\link{assoc}}, \code{\link{test_rho}}, \code{\link{test_beta}}
}
\examples{
  library(LPM)
  
  data(ExampleData)
  
  bLPMfit_first3 <- bLPM(ExampleData$data[1:3], X = ExampleData$X, coreNum = 2)
  bLPMfit <- bLPM_add(ExampleData$data[1:3], ExampleData$data[4], X = ExampleData$X, bLPMfit_first3, coreNum = 2)
  LPMfit <- LPM(bLPMfit)
  
  posterior1 <- post(ExampleData$data[1], X = ExampleData$X, id = 1, LPMfit)
  posterior13 <- post(ExampleData$data[c(1, 3)], X = ExampleData$X, id = c(1, 3), LPMfit)
  assoc.SNP <- assoc(posterior1, FDRset = 0.1, fdrControl = "global")
  
  p_value_test_rho <- test_rho(bLPMfit)
  result_test_beta <- test_beta(ExampleData$data, X = ExampleData$X, id = 1, LPMfit)
}
