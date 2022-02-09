\name{rFactorCopulaREMADA}
\Rdversion{1.1}
\alias{rFactorCopulaREMADA}
\alias{rFactorCopulaREMADA.norm}
\alias{rFactorCopulaREMADA.beta}
\title{
Simulation from  1-factor copula mixed models for joint meta-analysis of \eqn{T} diagnostic tests}

\description{
Simulation from  1-factor copula mixed models for joint meta-analysis of \eqn{T} diagnostic tests
}

\usage{
rFactorCopulaREMADA.norm(N,p,si,taus,qcond1,tau2par1,qcond2,tau2par2)
rFactorCopulaREMADA.beta(N,p,g,taus,qcond1,tau2par1,qcond2,tau2par2)      }

\arguments{
\item{N}{number of studies}
\item{p}{vector  of sensitivities and specificities}
\item{si}{vector of variabilities; normal margins}
\item{g}{vector of variabilities; beta margins}
\item{taus}{Kendall's tau values}
\item{qcond1}{function for the inverse conditional copula cdfs that link the factor with the latent sensitivities}
\item{tau2par1}{function for maping Kendall's tau  to copula parameter at the copulas that link the factor with the latent sensitivities}
\item{qcond2}{function for the inverse conditional copula cdfs that link the factor with the latent specificities}
\item{tau2par2}{function for maping Kendall's tau  to copula parameter at the copulas that link the factor with the latent specificities}
}

\value{A list with the simulated data in matrices with \eqn{T} columns and \eqn{N} rows. \describe{
    \item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
}
}
\references{
Nikoloulopoulos, A.K. (2020) An one-factor copula mixed model for joint meta-analysis of multiple diagnostic tests.  \emph{ArXiv e-prints}. arXiv:2006.09278. \url{https://arxiv.org/abs/2006.09278}.
}
\examples{
N=50

qcond1=qcondcln
tau2par1=tau2par.cln
qcond2=qcondcln270
tau2par2=tau2par.cln270

p=c(0.8,0.7,0.8,0.7,0.8,0.7)
mu=log(p/(1-p))
si=rep(1,6)
taus=c(0.6,0.7,0.5,-0.3,-0.4,-0.2)

out=rFactorCopulaREMADA.norm(N,p,si,taus,qcond1,tau2par1,qcond2,tau2par2)
  
TP=out$TP
FN=out$FN
TN=out$TN
FP=out$FP
}
\keyword{simulation}

