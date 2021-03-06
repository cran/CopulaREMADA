\name{dvinesim}
\Rdversion{1.1}
\alias{dvinesim}
\alias{dtrvinesim}
\title{
Simulation from a (truncated) quadrivariate D-vine copula
}
\description{
Simulation from a (truncated) quadrivariate D-vine copula. Lower-order trees (if any) are composed with BVN copulas. 
}
\usage{
dvinesim(nsim,param,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)
dtrvinesim(nsim,trparam,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)

}
\arguments{
\item{nsim}{sample size}
\item{param}{Kendall's tau values for margins (1,2), (2,3), (3,4), (1,3|2), (2,4|3), (1,4|23)}
\item{trparam}{Kendall's tau values for margins (1,2), (2,3), (3,4)}
\item{qcond1}{function for the inverse conditional copula cdf at the (1,2) and (3,4) bivariate margins}
\item{pcond1}{function for the conditional copula cdf at the (1,2) and  and (3,4) bivariate margins}
\item{tau2par1}{function for maping Kendall's tau at the (1,2) and (3,4) bivariate margins to copula parameter}
\item{qcond2}{function for the inverse conditional copula cdf at the (2,3) bivariate margin}
\item{pcond2}{function for the conditional copula cdf at the (2,3) bivariate margin}
\item{tau2par2}{function for maping Kendall's tau at the (2,3) bivariate margin to copula parameter}
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

Joe H (2014) \emph{Dependence Modeling with Copulas}. Chapman & Hall/CRC.

Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. Software for book: \emph{Dependence Modeling with Copulas}, Chapman & Hall/CRC, 2014.

Nikoloulopoulos, A.K. (2018) A D-vine copula mixed model for joint meta-analysis and comparison of  diagnostic tests.  \emph{Statistical Methods in Medical Research}, in press. \doi{10.1177/0962280218796685}. 

Nikoloulopoulos, A.K. (2018) A multinomial quadrivariate D-vine copula mixed model for diagnostic studies meta-analysis accounting for non-evaluable subjects.  \emph{ArXiv e-prints}, arXiv:1812.05915. \url{https://arxiv.org/abs/1812.05915}.
}
\seealso{
\code{\link{qcondcop}}
\code{\link{dcop}}
\code{\link{rcop}}
}


\keyword{copula}
\keyword{simulation}
