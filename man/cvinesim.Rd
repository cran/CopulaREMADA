\name{cvinesim}
\Rdversion{1.1}
\alias{cvinesim}
\title{
Simulation from a trivariate C-vine copula
}
\description{
Simulation from a trivariate C-vine copula
}
\usage{
cvinesim(N,param,qcondcop12,qcondcop13,qcondcop23,
                  tau2par12,tau2par13,tau2par23)
}
\arguments{
\item{N}{sample size}
\item{param}{Kendall's tau values for margins (1,2), (1,3), (23|1)}
\item{qcondcop12}{function for the inverse of conditional copula cdf at the (1,2) bivariate margin}
\item{qcondcop13}{function for the inverse of conditional copula cdf at the (1,3) bivariate margin}
\item{qcondcop23}{function for the inverse of conditional copula cdf at the (2,3|1) bivariate margin}
\item{tau2par12}{function for maping Kendall's tau at the (1,2) bivariate margin to copula parameter}
\item{tau2par13}{function for maping Kendall's tau at the (1,3) bivariate margin to copula parameter}
\item{tau2par23}{function for maping Kendall's tau at the (2,3|1) bivariate margin to the conditional copula parameter}
}
\details{
Choices are 'cop' in rcop are bvn, frk, cln, cln90 (rotated by 90 degrees cln), cln180 (rotated by 180 degrees cln), cln270 (rotated by 270 degrees cln).

See help page for \code{\link{dcop}} for the abbreviations of the copula names.

}
\value{
Nx3 matrix with values in (0,1)
}
\references{
Joe H (2011) Dependence comparisons of vine copulae with four or more variables. In: Kurowicka D, Joe H, editors. \emph{Dependence Modeling: Handbook on Vine Copulae}. Singapore: World Scientific; 2011. p. 139--164

Joe H (2014) \emph{Dependence Modeling with Copulas}. Chapman & Hall/CRC.

Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. Software for book: \emph{Dependence Modeling with Copulas}, Chapman & Hall/CRC, 2014.
}
\seealso{
\code{\link{qcondcop}}
\code{\link{dcop}}
\code{\link{rcop}}
}


\keyword{copula}
\keyword{simulation}
