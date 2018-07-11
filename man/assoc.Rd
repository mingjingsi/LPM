\name{assoc}
\alias{assoc}
\title{
  Statistical inference for association mapping
}
\description{
Statistical inference for association mapping.
}
\usage{
assoc(post, FDRset = 0.1, fdrControl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{post}{
     The posterior for association mapping.
}
\item{FDRset}{
     FDR level. Default is 0.1.
}
  \item{fdrControl}{
     Method to control FDR. Possible values are "global" (global FDR control) and "local" (local FDR control).
}
}
\details{
  \code{assoc} makes statistical inference for association mapping.
}
\value{
The binary values indicating association of SNPs, where one indicates association and zero otherwise. The format is the same with \code{post}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jingsi Ming
}
\seealso{
 \code{\link{post}}
}
\examples{
assoc.SNP <- assoc(posterior1, FDRset = 0.1, fdrControl = "global")
}