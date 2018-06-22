\name{dvinesim}
\Rdversion{1.1}
\alias{dvinesim}
\title{
Simulation from a quadrivariate D-vine copula
}
\description{
Simulation from a quadrivariate D-vine copula
}
\usage{
dvinesim(nsim,param,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)

}
\arguments{
\item{nsim}{sample size}
\item{param}{Kendall's tau values for margins (1,2), (2,3), (3,4), (1,3|2), (2,4|3), (1,4|23)}
\item{qcond1}{function for the inverse conditional copula cdf at the (1,2) bivariate margin}
\item{pcond1}{function for the conditional copula cdf at the (1,2) bivariate margin}
\item{tau2par1}{function for maping Kendall's tau at the (1,2) bivariate margin to copula parameter}
\item{qcond2}{function for the inverse conditional copula cdf at the (3,4) bivariate margin}
\item{pcond2}{function for the conditional copula cdf at the (3,4) bivariate margin}
\item{tau2par2}{function for maping Kendall's tau at the (3,4) bivariate margin to copula parameter}
}
\details{
Choices are 'cop' in rcop are bvn, frk, cln, cln90 (rotated by 90 degrees cln), cln180 (rotated by 180 degrees cln), cln270 (rotated by 270 degrees cln).

See help page for \code{\link{dcop}} for the abbreviations of the copula names.

}
\value{
Nx4 matrix with values in (0,1)
}
\references{
Joe H (2011) Dependence comparisons of vine copulae with four or more variables. In: Kurowicka D, Joe H, editors. \emph{Dependence Modeling: Handbook on Vine Copulae}. Singapore: World Scientific; 2011. p. 139--164

Joe H (2014) \emph{Dependence Modeling with Copulas}. Chapman & Hall

Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. \url{http://copula.stat.ubc.ca/}. 

Nikoloulopoulos, A.K.  (2018b) A D-vine copula mixed model for joint meta-analysis and comparison of  diagnostic tests.  \emph{ArXiv e-prints}, arXiv:1805.09674. \url{https://arxiv.org/abs/1805.09674}. 
}
\seealso{
\code{\link{qcondcop}}
\code{\link{dcop}}
\code{\link{rcop}}
}


\keyword{copula}
\keyword{simulation}
