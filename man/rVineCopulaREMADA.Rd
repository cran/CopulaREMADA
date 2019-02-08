\name{rVineCopulaREMADA}
\Rdversion{1.1}
\alias{rVineCopulaREMADA}
\alias{rVineCopulaREMADA.norm}
\alias{rVineCopulaREMADA.beta}
\title{
Simulation from trivariate vine copula mixed models for diagnostic test accuaracy studies accounting for disease prevalence and non-evaluable results
}

\description{
Simulation from trivariate vine copula mixed models for diagnostic test accuaracy studies accounting for disease prevalence and non-evaluable results
}

\usage{
rVineCopulaREMADA.beta(N,p,g,taus,omega1,omega0,qcondcop12,qcondcop13,
qcondcop23,tau2par12,tau2par13,tau2par23)
rVineCopulaREMADA.norm(N,p,si,taus,omega1,omega0,qcondcop12,qcondcop13,
qcondcop23,tau2par12,tau2par13,tau2par23)
}
\arguments{
\item{N}{sample size}
\item{p}{Vector \eqn{(\pi_1,\pi_2,\pi_3)} of sensitivity/specificity/prevalence}
\item{si}{Vector \eqn{(\sigma_1,\sigma_2,\sigma_3)} of variability; normal margins}
\item{g}{Vector \eqn{(\gamma_1,\gamma_2,\gamma_3)} of variability; beta margins}
\item{taus}{Kendall's tau values}
\item{omega1}{the probability for non-evaluable positives}
\item{omega0}{the probability for non-evaluable negatives}
\item{qcondcop12}{function for the inverse of conditional copula cdf at the (1,2) bivariate margin}
\item{qcondcop13}{function for the inverse of conditional copula cdf at the (1,3) bivariate margin}
\item{qcondcop23}{function for the inverse of conditional copula cdf at the (2,3|1) bivariate margin}
\item{tau2par12}{function for maping Kendall's tau at the (1,2) bivariate margin to copula parameter}
\item{tau2par13}{function for maping Kendall's tau at the (1,3) bivariate margin to copula parameter}
\item{tau2par23}{function for maping Kendall's tau at the (2,3|1) bivariate margin to the conditional copula parameter}
}

\value{Simuated data with 6 columns and \eqn{N} rows. \describe{
    \item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
\item{NEP}{the number of non-evaluable positives}
\item{NEN}{the number of non-evaluable negatives}
}
}
\references{
Nikoloulopoulos, A.K.  (2017) A vine copula mixed effect model for trivariate meta-analysis of diagnostic test accuracy studies accounting for disease prevalence. \emph{Statistical Methods in Medical Research}, \bold{26}, 2270--2286. \doi{10.1177/0962280215596769}.

Nikoloulopoulos, A.K.  (2018) A vine copula mixed  model for trivariate meta-analysis of diagnostic studies accounting for disease prevalence and non-evaluable subjects.  \emph{ArXiv e-prints}, arXiv:1812.03685. \url{https://arxiv.org/abs/1812.03685}.
}
\seealso{
\code{\link{rCopulaREMADA}} \code{\link{rcop}}  \code{\link{cvinesim}}
}
\examples{
p=c(0.8,0.7,0.4)
g=c(0.1,0.1,0.05)
taus=c(-0.5,-0.3,-0.0001)
qcondcop12=qcondcop23=qcondcop13=qcondcln90
tau2par12=tau2par23=tau2par13=tau2par.cln90
# in the absence of non-evaluable results
omega1=0
omega0=0
rVineCopulaREMADA.beta(50,p,g,taus,omega1,omega0,
qcondcop12,qcondcop13,qcondcop23,tau2par12,
tau2par13,tau2par23)
# in the presence of non-evaluable results
omega1=0.1
omega0=0.2
rVineCopulaREMADA.beta(50,p,g,taus,omega1,omega0,
qcondcop12,qcondcop13,qcondcop23,tau2par12,
tau2par13,tau2par23)
}
\keyword{simulation}

