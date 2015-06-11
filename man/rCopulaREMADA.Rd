\name{rCopulaREMADA}
\Rdversion{1.1}
\alias{rCopulaREMADA}
\alias{rCopulaREMADA.norm}
\alias{rCopulaREMADA.beta}
\title{
Simulation from copula mixed models for diagnostic test accuaracy studies
}
\description{
To simulate the data we have used the following steps:

1. Simulate the study size \eqn{n} from a shifted gamma distribution with parameters \eqn{\alpha=1.2,\beta=0.01,lag=30} and round off to the nearest integer. 

2. Simulate \eqn{(u_1,u_2)} from a parametric family of copulas 'cop'.  

3. Convert to beta realizations  or normal realizations. 

4. Draw the number of diseased \eqn{n_1} from a \eqn{B(n,0.43)} distribution.

5. Set  \eqn{n_2=n-n_1, y_j=n_jx_j} and then round \eqn{y_j} for \eqn{j=1,2}. 

 
}
\usage{
rCopulaREMADA.norm(N,p,si,tau,rcop,tau2par)
rCopulaREMADA.beta(N,p,g,tau,rcop,tau2par) 
}
\arguments{
\item{N}{sample size}
\item{p}{Pair \eqn{(\pi_1,\pi_2)} of sensitivity/specificity}
\item{si}{Pair \eqn{(\sigma_1,\sigma_2)} of variability; normal margins}
\item{g}{Pair \eqn{(\gamma_1,\gamma_2)} of variability; beta margins}
\item{tau}{Kendall's tau value}
\item{rcop}{function for copula generation}
\item{tau2par}{function for mapping from Kendall's tau to copula parameter }
}
\value{A list containing the following simulated components:
\item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
}
\references{
Nikoloulopoulos, A.K. (2015a) A mixed effect model for bivariate meta-analysis of diagnostic test accuracy studies using  a copula representation of the random effects distribution, \emph{ArXiv e-prints}.
}
\seealso{
\code{\link{CopulaREMADA}} \code{\link{rcop}}
}
\examples{
nq=15
gl=gauss.quad.prob(nq,"uniform")
mgrid<- meshgrid(gl$n,gl$n)

N=20
tau=-0.5
p=c(0.7,0.9)
g=c(0.2,0.1)
simDat=rCopulaREMADA.beta(N,p,g,tau,rcln270,tau2par.cln270)
TP=simDat$TP
TN=simDat$TN
FP=simDat$FP
FN=simDat$FN
c270est.b=CopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid,qcondcln270,tau2par.cln270)

si=c(2,1)
tau=0.5
simDat=rCopulaREMADA.norm(N,p,si,tau,rcln,tau2par.cln)
TP=simDat$TP
TN=simDat$TN
FP=simDat$FP
FN=simDat$FN
cest.n=CopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid,qcondcln,tau2par.cln)
}
\keyword{simulation}

