\name{dcop}
\Rdversion{1.1}
\alias{dcop}
\alias{dbvn}
\alias{dfrk}
\alias{dcln}
\alias{dcln90}
\alias{dcln180}
\alias{dcln270}
\title{
Bivariate copula densities
}
\description{
Bivariate copula densities for parametric families.
}
\usage{
dbvn(u,v,cpar)
dfrk(u,v,cpar)
dcln(u,v,cpar)
dcln90(u,v,cpar)
dcln270(u,v,cpar)
}
\arguments{
\item{u}{value in interval 0,1; could be a vector}
\item{v}{value in interval 0,1; could be a vector}
\item{cpar}{copula parameter: scalar.}
}
\details{
Choices are 'cop' in dcop are bvn, frk, cln, cln90 (rotated by 90 degrees cln), 
cln180 (rotated by 180 degrees cln), cln270 (rotated by 270 degrees cln).

The copula names are abbreviations for:

bvn = bivariate normal or Gaussian

frk = Frank

cln = Clayton or Mardia-Takahasi-Cook-Johnson

}
\value{
pdf value(s).
}
\references{
Joe H (1997) \emph{Multivariate Models and Dependence Concepts}. Chapman & Hall

Joe H (2014) \emph{Dependence Modeling with Copulas}. Chapman & Hall

Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. \url{http://copula.stat.ubc.ca/}. 
}
\seealso{
\code{\link{qcondcop}}
\code{\link{rcop}}
}


\keyword{copula}
\keyword{distribution}
