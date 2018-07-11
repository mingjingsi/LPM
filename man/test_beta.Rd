\name{test_beta}
\alias{test_beta}
\title{
  Hypothesis testing of annotation enrichment
}
\description{
Hypothesis testing of annotation enrichment.
}
\usage{
test_beta(data, X, id, LPMfit)
test_beta_louise(data, X, id, LPMfit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     A list of dataframes which contains SNP IDs and p-values for multiple GWASs, the length is the number of GWASs.
}
  \item{X}{
     Design matrix of functional annotations without intercept, where row and column correspond to SNP and annotation, respectively. Default is \code{NULL}.
}
  \item{id}{
     ID of the target GWAS, the length is 1. For example, \code{id = 1} means we consider the first GWAS in \code{data}.
}
  \item{LPMfit}{
     LPM model fit.
}
}
\details{
  \code{test_beta} implements the hypothesis testing of annotation enrichment.
  \code{test_beta_louise} implements the hypothesis testing of annotation enrichment using Louise method.
}
\value{
Returns 
  \item{p_value}{
     p-value of the hypothesis testing of annotation enrichment.
}
  \item{se}{
     Standard error of the estimates \code{beta}.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming
}
\seealso{
 \code{\link{LPM}}, \code{\link{test_beta}}
}
\examples{
data(ExampleData)

result_test_beta <- test_beta(ExampleData$data, X = ExampleData$X, id = 1, LPMfit)

}