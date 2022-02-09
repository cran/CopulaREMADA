\name{rmultinomVineCopulaREMADA}
\Rdversion{1.1}
\alias{rmultinomVineCopulaREMADA}
\alias{rmultinomVineCopulaREMADA.norm}
\alias{rmultinomVineCopulaREMADA.beta}
\title{
Simulation from  multinomial quadrivariate (truncated) D-vine copula mixed models for diagnostic test accurracy studies accounting for non-evaluable outcomes
}

\description{
Simulation from   multinomial quadrivariate (truncated) D-vine copula mixed models for diagnostic test accurracy studies accounting for non-evaluable outcomes
}

\usage{
rmultinomVineCopulaREMADA.norm(N,p,si,taus,qcond1,
                                    pcond1,tau2par1,qcond2,
                                    pcond2,tau2par2)
rmultinomVineCopulaREMADA.beta(N,p,g,taus,qcond1,
                                    pcond1,tau2par1,qcond2,
                                    pcond2,tau2par2)                                    
}
\arguments{
\item{N}{sample size}
\item{p}{Vector \eqn{(\pi_1,\pi_2,\pi_3)} of sensitivity/specificity/prevalence}
\item{si}{Vector \eqn{(\sigma_1,\sigma_2,\sigma_3)} of variability; normal margins}
\item{g}{Vector \eqn{(\gamma_1,\gamma_2,\gamma_3)} of variability; beta margins}
\item{taus}{Kendall's tau values}
\item{qcond1}{function for the inverse conditional copula cdf at the (1,2) and (3,4) bivariate margin}
\item{pcond1}{function for the conditional copula cdf at the (1,2) and (3,4) bivariate margin}
\item{tau2par1}{function for maping Kendall's tau at the (1,2) and (3,4) bivariate margin to copula parameter}
\item{qcond2}{function for the inverse conditional copula cdf at the (2,3) bivariate margin}
\item{pcond2}{function for the conditional copula cdf at the (2,3) bivariate margin}
\item{tau2par2}{function for maping Kendall's tau at the (2,3) bivariate margin to copula parameter}
}

\value{Simulated data with 6 columns and \eqn{N} rows. \describe{
    \item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
\item{NEP}{the number of non-evaluable positives}
\item{NEN}{the number of non-evaluable negatives}
}
}
\references{
Nikoloulopoulos, A.K. (2020) A multinomial quadrivariate D-vine copula mixed model for diagnostic studies meta-analysis in the presence of non-evaluable subjects.  \emph{Statistical Methods in Medical Research}, 29 (10), 2988--3005. \doi{10.1177/0962280220913898}.

}
\seealso{
\code{\link{dvinesim}}}

\examples{
N=30
p=c(0.898745016,0.766105342,0.059168715,0.109217888)
g=c(0.090270947,0.079469009,0.367463579,0.154976269)
taus=c(	0.82050793,-0.51867629,0.26457961)

qcond1=qcondcln180
pcond1=pcondcln180
tau2par1=tau2par.cln180

qcond2=qcondcln90
pcond2=pcondcln90
tau2par2=tau2par.cln90

out=rmultinomVineCopulaREMADA.beta(N,p,g,taus,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)
  
TP=out[,1]
NEP=out[,2]
FN=out[,3]
TN=out[,4]
NEN=out[,5]
FP=out[,6]
}
\keyword{simulation}

