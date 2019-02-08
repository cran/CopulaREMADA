\name{rcop}
\Rdversion{1.1}
\alias{rcop}
\alias{rbvn}
\alias{rfrk}
\alias{rcln}
\alias{rcln90}
\alias{rcln180}
\alias{rcln270}
\title{
Simulation from parametric bivariate copula families
}
\description{
Simulation from parametric bivariate copula families
}
\usage{
rbvn(N,cpar)
rfrk(N,cpar)
rcln(N,cpar)
rcln90(N,cpar)
rcln270(N,cpar) 
}
\arguments{
\item{N}{sample size}
\item{cpar}{copula parameter: scalar}
}
\details{
Choices are 'cop' in rcop are bvn, frk, cln, cln90 (rotated by 90 degrees cln), cln180 (rotated by 180 degrees cln), cln270 (rotated by 270 degrees cln).

See help page for \code{\link{dcop}} for the abbreviations of the copula names.

}
\value{
nx2 matrix with values in (0,1)
}
\references{
Joe H (1997) \emph{Multivariate Models and Dependence Concepts}. Chapman & Hall

Joe H (2014) \emph{Dependence Modeling with Copulas}. Chapman & Hall/CRC.

Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. Software for book: \emph{Dependence Modeling with Copulas}, Chapman & Hall/CRC, 2014. 
}
\seealso{
\code{\link{qcondcop}}
\code{\link{dcop}}
}


\keyword{copula}
\keyword{simulation}
