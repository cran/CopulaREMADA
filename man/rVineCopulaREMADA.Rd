\name{rVineCopulaREMADA}
\Rdversion{1.1}
\alias{rVineCopulaREMADA}
\alias{rVineCopulaREMADA.norm}
\alias{rVineCopulaREMADA.beta}
\title{
Simulation from vine copula mixed models for diagnostic test accuaracy studies accounting for disease prevalence
}
\description{
To simulate the data we have used the following steps:

1. Simulate the study size \eqn{n} from a shifted gamma distribution with parameters \eqn{\alpha=1.2,\beta=0.01,lag=30} and round off to the nearest integer. 

2. Simulate \eqn{(u_1,u_2,u_3)} from a trivariate C-vine copula. See also \code{\link{cvinesim}}.  

3. Convert to beta realizations  or normal realizations. 

4. Set the number of diseased \eqn{n_1=nx_1}.

5. Set  \eqn{n_2=n-n_1, y_j=n_jx_j} and then round \eqn{y_j} for \eqn{j=1,2}. 

 
}
\usage{
rVineCopulaREMADA.beta(N,p,g,taus,qcondcop12,qcondcop13,qcondcop23,
             tau2par12,tau2par13,tau2par23)
rVineCopulaREMADA.norm(N,p,si,taus,qcondcop12,qcondcop13,qcondcop23,
             tau2par12,tau2par13,tau2par23)
}
\arguments{
\item{N}{sample size}
\item{p}{Vector \eqn{(\pi_1,\pi_2,\pi_3)} of sensitivity/specificity/prevalence}
\item{si}{Vector \eqn{(\sigma_1,\sigma_2,\sigma_3)} of variability; normal margins}
\item{g}{Vector \eqn{(\gamma_1,\gamma_2,\gamma_3)} of variability; beta margins}
\item{taus}{Kendall's tau values}
\item{qcondcop12}{function for the inverse of conditional copula cdf at the (1,2) bivariate margin}
\item{qcondcop13}{function for the inverse of conditional copula cdf at the (1,3) bivariate margin}
\item{qcondcop23}{function for the inverse of conditional copula cdf at the (2,3|1) bivariate margin}
\item{tau2par12}{function for maping Kendall's tau at the (1,2) bivariate margin to copula parameter}
\item{tau2par13}{function for maping Kendall's tau at the (1,3) bivariate margin to copula parameter}
\item{tau2par23}{function for maping Kendall's tau at the (2,3|1) bivariate margin to the conditional copula parameter}
}

\value{A list containing the following simulated components:
\item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
}
\references{
Nikoloulopoulos, A.K.  (2017) A vine copula mixed effect model for trivariate meta-analysis of diagnostic test accuracy studies accounting for disease prevalence. \emph{Statistical Methods in Medical Research}, \bold{26}, 2270--2286. \doi{10.1177/0962280215596769}.
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
rVineCopulaREMADA.beta(50,p,g,taus,qcondcop12,qcondcop13,qcondcop23,
             tau2par12,tau2par13,tau2par23)
}
\keyword{simulation}

