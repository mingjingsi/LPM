\name{single_LPM}
\alias{single_LPM}
\title{
  Fit single_LPM
}
\description{
Fit single_LPM.
}
\usage{
single_LPM(data, X = NULL, alpha = 0.1, pi1_ = 0.1, maxiter = 1e4, tol = 1e-8, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     A dataframe which contains SNP IDs and p-values for single GWASs.
}
  \item{X}{
     Design matrix of functional annotations without intercept, where row and column correspond to SNP and annotation, respectively. Default is \code{NULL}.
}
  \item{alpha}{
     Initial value of parameter \code{alpha} in the Beta distribution which is between 0 and 1 for GWAS. Default is 0.1.
}
  \item{pi1_}{
     Initial value of parameter \code{pi1} which is the proportion that eta equals 1 for GWAS. Default is 0.1.
}
  \item{maxiter}{
     Maximum number of EM iteration. Default is 1000.
}
  \item{tol}{
     Threshold to stop the EM iteration. Default is 1e-08.
}
  \item{verbose}{
     Whether or not to output the estimates of parameters in stage 1. Defalut is FALSE.
  }
  
}
\details{
  It requires to provide SNP IDs and p-values for single GWASs to \code{data}, while users can also provide annotation data to \code{X}. The number of SNPs for GWASs are required to be the same.
}
\value{
List of parameter estimation, posterior and the value of lower bound of Log-likelihood.
  \item{alpha}{
     Estimation of parameter \code{alpha}.
}
  \item{alpha_stage1}{
     Estimation of parameter \code{alpha} in the Stage 1.
}
  \item{pi1}{
     The posteriors that each SNP is associated with the GWAS.
}
  \item{pi1_stage1}{
     The posteriors that each SNP is associated with the GWAS in Stage 1.
}
  \item{beta}{
     Estimation of parameter \code{beta} for each GWAS in each pair.
}
  \item{LL}{
     The value of Log-likelihood.
}
  \item{L_stage1_List}{
     The value of Log-likelihood in Stage 1.
}
  \item{L_stage2_List}{
     The value of Log-likelihood in Stage 2.
}


}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming
}
\seealso{
 \code{\link{LPM}}, \code{\link{bLPM}}
}
\examples{
data(ExampleData)

single_LPMfit <- single_LPM(data = ExampleData$data$P1, X = ExampleData$X)

}