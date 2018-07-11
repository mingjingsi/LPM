\name{bLPM_add}
\alias{bLPM_add}
\title{
  Fit bLPM if more GWASs are added
}
\description{
Fit bLPM if more GWASs are added.
}
\usage{
bLPM_add(data, data_add, X = NULL, bLPMfit, maxiter = 1e4, tol = 1e-8, coreNum = 1, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     A list of dataframes which contains SNP IDs and p-values for original GWASs, the length is the number of GWASs.
}
  \item{data}{
     A list of dataframes which contains SNP IDs and p-values for added GWASs, the length is the number of GWASs.
}
  \item{X}{
     Design matrix of functional annotations without intercept, where row and column correspond to SNP and annotation, respectively. Default is \code{NULL}.
}
  \item{bLPMfit}{
     bLPM fit for original GWASs.
}
  \item{maxiter}{
     Maximum number of EM iteration. Default is 1000.
}
  \item{tol}{
     Threshold to stop the EM iteration. Default is 1e-08.
}
  \item{coreNum}{
     Number of cores used. Default is 1.
  }
  \item{verbose}{
     Whether or not to output the estimates of parameters in stage 1 and stage2. Defalut is FALSE.
  }
  
}
\details{
  \code{bLPM_add} fits the bLPM when more GWASs are added.
}
\value{
List of parameter estimation, the value of lower bound of Log-likelihood. It has the same form as \code{bLPM}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming
}
\seealso{
 \code{\link{bLPM}}
}
\examples{
data(ExampleData)

bLPMfit <- bLPM(ExampleData$data[1:3], ExampleData$data[4], X = ExampleData$X, bLPMfit_first3, coreNum = 2)
}