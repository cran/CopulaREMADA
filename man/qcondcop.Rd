\name{qcondcop}
\Rdversion{1.1}
\alias{qcondcop}
\alias{qcondbvn}
\alias{qcondfrk}
\alias{qcondcln}
\alias{qcondcln90}
\alias{qcondcln180}
\alias{qcondcln270}
\title{
Bivariate copula conditional quantile functions
}
\description{
Bivariate copula conditional quantile functions
}
\usage{
qcondbvn(p,u,cpar)
qcondfrk(p,u,cpar)
qcondcln(p,u,cpar)
qcondcln90(p,u,cpar)
qcondcln270(p,u,cpar)
}
\arguments{
\item{u}{conditioning value in interval 0,1; could be a vector }
\item{p}{quantile in interval 0,1; could be a vector }
\item{cpar}{copula parameter: scalar. }
}
\details{
Choices appending 'cop'
are bvn, frk, cln, cln90 (rotated by 90 degrees cln), cln180 (rotated by 180 degrees cln), cln270 (rotated by 270 degrees cln).

See help page for \code{\link{dcop}} for the abbreviations of the copula names.
}
\value{inverse conditional cdf value(s)
}
\references{
Joe H (1997) \emph{Multivariate Models and Dependence Concepts}. Chapman & Hall

Joe H (2014) \emph{Dependence Modeling with Copulas}. Chapman & Hall/CRC.

Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. Software for book: \emph{Dependence Modeling with Copulas}, Chapman & Hall/CRC, 2014. 
}
\seealso{
\code{\link{dcop}}
\code{\link{rcop}}
}



\keyword{copula}
\keyword{distribution}
