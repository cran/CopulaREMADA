\name{pcondcop}
\Rdversion{1.1}
\alias{pcondcop}
\alias{pcondbvn}
\alias{pcondfrk}
\alias{pcondcln}
\alias{pcondcln90}
\alias{pcondcln180}
\alias{pcondcln270}
\title{
Bivariate copula conditional distribution functions
}
\description{
Bivariate copula conditional distribution functions
}
\usage{
pcondbvn(v,u,cpar)
pcondfrk(v,u,cpar)
pcondcln(v,u,cpar)
pcondcln90(v,u,cpar)
pcondcln270(v,u,cpar)
}
\arguments{
\item{v}{conditioning value in interval 0,1; could be a vector }
\item{u}{value in interval 0,1; could be a vector }
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

Joe H (2014) \emph{Dependence Modeling with Copulas}. Chapman & Hall

Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. \url{http://copula.stat.ubc.ca/}. 
}
\seealso{
\code{\link{dcop}}
\code{\link{rcop}}
}



\keyword{copula}
\keyword{distribution}
