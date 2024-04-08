\name{dvine6dsim}
\Rdversion{1.1}
\alias{dvine6dsim}
\title{
Simulation from a six-variate 1-truncated  D-vine copula
}
\description{
Simulation from a six-variate 1-truncated  D-vine copula.
}
\usage{
dvine6dsim(nsim,tau,qcond,tau2par)
}
\arguments{
\item{nsim}{sample size}
\item{tau}{Kendall's tau values}
\item{qcond}{function for the inverse conditional copula cdf}
\item{tau2par}{function for maping Kendall's taus to copula parameters}
}
\details{
Choices are 'cop' in rcop are bvn, frk, cln, cln90 (rotated by 90 degrees cln), cln180 (rotated by 180 degrees cln), cln270 (rotated by 270 degrees cln).

See help page for \code{\link{qcondcop}} for the abbreviations of the copula names.

}
\value{
Nx6 matrix with values in (0,1)
}
\references{
Joe H (2011) Dependence comparisons of vine copulae with four or more variables. In: Kurowicka D, Joe H, editors. \emph{Dependence Modeling: Handbook on Vine Copulae}. Singapore: World Scientific; 2011. p. 139--164

Joe H (2014) \emph{Dependence Modeling with Copulas}. Chapman & Hall/CRC.

Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. Software for book: \emph{Dependence Modeling with Copulas}, Chapman & Hall/CRC, 2014.

Nikoloulopoulos, A.K. (2024) Joint meta-analysis  of two diagnostic tests accounting for within and between studies dependence. Submitted.

}
\seealso{
\code{\link{qcondcop}}
\code{\link{rcop}}
}


\keyword{copula}
\keyword{simulation}
