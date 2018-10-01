\name{post}
\alias{post}
\title{
  Compute posterior for association mapping
}
\description{
Compute posterior for association mapping.
}
\usage{
post(data, X = NULL, id, LPMfit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     A list of dataframes which contains SNP IDs and p-values for the target GWASs, the length is the number of GWASs considered.
}
  \item{X}{
     Design matrix of functional annotations without intercept, where row and column correspond to SNP and annotation, respectively. Default is \code{NULL}.
}
  \item{id}{
     ID of the GWASs considered, the length can be 1, 2 and 3. For example, \code{id = c(1, 3)} means we consider the first and the third GWAS in \code{data}.
}
  \item{LPMfit}{
     LPM model fit.
}
}
\details{
  \code{post} compute posterior for association mapping.
}
\value{
If the length of \code{id} is 1, the value includes
  \item{posterior}{
     The posteriors that each SNP is associated with the target GWAS.
}
If the length of \code{id} is 2, the value includes
  \item{post.joint}{
     The posteriors that each SNP is associated with both the two target GWASs.
}
\item{post.marginal1}{
     The posteriors that each SNP is associated with only the first target GWAS.
}
\item{post.marginal2}{
     The posteriors that each SNP is associated with only the second target GWAS.
}
If the length of \code{id} is 2, the value includes
  \item{post.joint}{
     The posteriors that each SNP is associated with all the three target GWASs.
}
\item{post.marginal1}{
     The posteriors that each SNP is associated with only the first target GWAS.
}
\item{post.marginal2}{
     The posteriors that each SNP is associated with only the second target GWAS.
}
\item{post.marginal3}{
     The posteriors that each SNP is associated with only the third target GWAS.
}
\item{post.marginal12}{
     The posteriors that each SNP is associated with both the first and the second target GWAS.
}
\item{post.marginal13}{
     The posteriors that each SNP is associated with both the first and the third target GWAS.
}
\item{post.marginal23}{
     The posteriors that each SNP is associated with both the second and the third target GWAS.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming
}
\seealso{
 \code{\link{LPM}}, \code{\link{assoc}}
}
\examples{
data(ExampleData)

posterior1 <- post(ExampleData$data[1], X = ExampleData$X, id = 1, LPMfit)

posterior13 <- post(ExampleData$data[c(1, 3)], X = ExampleData$X, id = c(1, 3), LPMfit)
}