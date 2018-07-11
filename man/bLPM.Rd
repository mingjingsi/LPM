\name{bLPM}
\alias{bLPM}
\title{
  Fit bLPM
}
\description{
Fit bLPM.
}
\usage{
bLPM(data, X = NULL, alpha = NULL, pi1_ = NULL, maxiter = 1e4, tol = 1e-8, coreNum = 1, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     A list of dataframes which contains SNP IDs and p-values for multiple GWASs, the length is the number of GWASs.
}
  \item{X}{
     Design matrix of functional annotations without intercept, where row and column correspond to SNP and annotation, respectively. Default is \code{NULL}.
}
  \item{alpha}{
     Initial values of parameter \code{alpha} in the Beta distribution which is between 0 and 1 for GWASs, the length is the number of GWASs. Default is 0.1.
}
  \item{pi1_}{
     Initial value of parameter \code{pi1} which is the proportion that eta equals 1 for GWASs, the length is the number of GWASs. Default is 0.1.
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
  \code{bLPM} fits the bLPM. It requires to provide SNP IDs and p-values for multiple GWASs to \code{data}, while users can also provide annotation data to \code{X}. The number of SNPs for GWASs are not required to be the same.
}
\value{
List of parameter estimation, the value of lower bound of Log-likelihood.
  \item{alpha}{
     Estimation of parameter \code{alpha} for each GWAS in each pair.
}
  \item{alpha_stage1}{
     Estimation of parameter \code{alpha} in the Stage 1 for each GWAS in each pair.
}
  \item{alpha_stage2}{
     Estimation of parameter \code{alpha} in the Stage 2 for each GWAS in each pair.
}
  \item{pi1_stage1}{
     Estimation of parameter pi1 in Stage 1 for each GWAS in each pair, which measures the proportion that eta equals 1.
}
  \item{beta0}{
     Estimation of parameter \code{beta} for each GWAS in each pair when annotation is not provided.
}
  \item{beta0_stage2}{
     Estimation of parameter \code{beta} in Stage 2 for each GWAS in each pair when annotation is not provided.
}
  \item{beta}{
     Estimation of parameter \code{beta} for each GWAS in each pair.
}
  \item{beta_stage2}{
     Estimation of parameter \code{beta} in Stage 2 for each GWAS in each pair.
}
  \item{rho}{
     Estimation of parameter \code{rho} for each pair of GWASs.
}
  \item{R}{
     Reform \code{rho} to the correlation matrix \code{R}.
}
  \item{L_stage1_List}{
     The value of Log-likelihood in Stage 1.
}
  \item{L_stage2_List}{
     The value of Log-likelihood in Stage 2.
}
  \item{L_stage3_List}{
     The value of Log-likelihood in Stage 3.
}


}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming
}
\seealso{
 \code{\link{LPM}}, \code{\link{test_rho}}
}
\examples{
data(ExampleData)

bLPMfit_first3 <- bLPM(ExampleData$data[1:3], X = ExampleData$X, coreNum = 2)

}