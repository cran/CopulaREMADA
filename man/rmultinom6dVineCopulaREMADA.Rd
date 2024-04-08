\name{rmultinom6dVineCopulaREMADA}
\Rdversion{1.1}
\alias{rmultinom6dVineCopulaREMADA}
\alias{rmultinom6dVineCopulaREMADA.norm}
\alias{rmultinom6dVineCopulaREMADA.beta}
\title{
Simulation from  multinomial six-variate 1-truncated D-vine copula mixed models for meta-analysis of two diagnostic tests accounting for within and between studies dependence
}

\description{
Simulation from   multinomial six-variate 1-truncated D-vine copula mixed models for meta-analysis of two diagnostic tests accounting for within and between studies dependence
}

\usage{rmultinom6dVineCopulaREMADA.norm(N,p,si,taus,qcond,tau2par)
rmultinom6dVineCopulaREMADA.beta(N,p,g,taus,qcond,tau2par)}
\arguments{
\item{N}{sample size}
\item{p}{Vector \eqn{(\pi_{101},\pi_{011},\pi_{111},\pi_{100},\pi_{010},\pi_{110})} of the meta-analytic parameters of interest for each combination of test results in diseased and non-diseased participants in a \eqn{4\times 2} table}
\item{si}{Vector \eqn{(\sigma_{101},\sigma_{011},\sigma_{111},\sigma_{100},\sigma_{010},\sigma_{110})} of variability parameters; normal margins}
\item{g}{Vector \eqn{(\gamma_{101},\gamma_{011},\gamma_{111},\gamma_{100},\gamma_{010},\gamma_{110})} of variability parameters; beta margins}
\item{taus}{Kendall's tau values}
\item{qcond}{function for the inverse conditional copula cdf}
\item{tau2par}{function for maping Kendall's taus to copula parameters}
}

\value{Simulated data with 8 columns and \eqn{N} rows. \describe{
    \item{y001}{the number of the test results in the diseased where the test 1 outcome is negative and  the test 2 outcome is negative}
\item{y011}{the number of the test results in the diseased where the test 1 outcome is negative and  the test 2 outcome is positive}
\item{y101}{the number of the test results in the diseased where the test 1 outcome is positive and  the test 2 outcome is negative}
\item{y111}{the number of the test results in the diseased where the test 1 outcome is positive and  the test 2 outcome is positive}
\item{y000}{the number of the test results in the non-diseased where the test 1 outcome is negative and  the test 2 outcome is negative}
\item{y010}{the number of the test results in the non-diseased where the test 1 outcome is negative and  the test 2 outcome is positive}
\item{y100}{the number of the test results in the non-diseased where the test 1 outcome is positive and  the test 2 outcome is negative}
\item{y110}{the number of the test results in the non-diseased where the test 1 outcome is positive and  the test 2 outcome is positive}
}
}
\references{
Nikoloulopoulos, A.K. (2024) Joint meta-analysis  of two diagnostic tests accounting for within and between studies dependence. Submitted.
}
\seealso{
\code{\link{dvine6dsim}}}

\examples{
N=11
p=c(0.03667409,  0.09299767,  0.29450436,  0.01733081,  0.04923809,  0.02984361)
si=c(1.69868880, 0.54292079,  0.58489574,  0.92918177,  0.48998484,  0.57004098)
taus=c(-0.52475006,  0.55768873, 0.18454559,  0.02233204,  0.57570506)


tau2par=tau2par.bvn
qcond=qcondbvn

out=rmultinom6dVineCopulaREMADA.norm(N,p,si,taus,qcond,tau2par)
 
  
y101=out[,1]
y011=out[,2]
y111=out[,3]
y001=out[,4]
y100=out[,5]
y010=out[,6]
y110=out[,7]
y000=out[,8] 
}
\keyword{simulation}

