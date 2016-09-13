\name{hybridCopulaREMADA}
\Rdversion{1.1}
\alias{hybridCopulaREMADA}
\alias{hybridCopulaREMADA.norm}
\alias{hybridCopulaREMADA.beta}


\title{
Maximum likelhood estimation  for hybrid copula mixed models for combining case-control and cohort studies in meta-analysis of diagnostic tests
}
\description{
The estimated parameters can be obtained by 
using a quasi-Newton method applied to the logarithm of the joint likelihood.  
This numerical  method requires only the objective
function, i.e.,  the logarithm of the joint likelihood, while the gradients
are computed numerically and the Hessian matrix of the second
order derivatives is updated in each iteration. The standard errors (SE) of the ML estimates can be also obtained via the gradients and the Hessian computed numerically during the maximization process. 
}
\usage{
hybridCopulaREMADA.norm(TP,FN,FP,TN,type,perm,gl,mgrid1,mgrid2,
               qcondcop12,qcondcop13,
               tau2par12,tau2par13,qcond,tau2par)
hybridCopulaREMADA.beta(TP,FN,FP,TN,type,perm,gl,mgrid1,mgrid2,
               qcondcop12,qcondcop13,
               tau2par12,tau2par13,qcond,tau2par) 
}
\arguments{
\item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
\item{type}{a scalar indicating the study type: 
1: Cohort study;

2: Case-control study.}
\item{perm}{a scalar indicating the selected permutation of indices: 

1: Pilot variable is the number of TP. The bivariate margins are 12, 13, 23|1; 

2: Pilot variable is the number of TN. The bivariate margins are 23, 12, 13|2; 

3: Pilot variable is the TP+FN. The bivariate margins are 13, 23, 12|3; 

1:TP, 2:TN, 3:TP+FN}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid1}{
a list containing three-dimensional arrays. For more details see \code{\link{meshgrid}}}
\item{mgrid2}{
a list containing two matrices with the rows of the output matrix x are copies of the vector \code{gl$nodes}; columns of the output matrix y are copies of the vector \code{gl$nodes}. For more details see also \code{\link{meshgrid}}}
\item{qcondcop12}{function for conditional copula cdf at the (1,2) bivariate margin of the vine}
\item{qcondcop13}{function for conditional copula cdf at the (1,3) bivariate margin of the vine}
\item{tau2par12}{function for maping Kendall's tau at the (1,2) bivariate margin of the vine to copula parameter}
\item{tau2par13}{function for maping Kendall's tau at the (1,3) bivariate margin of the vine to copula parameter}
\item{qcond}{function for the bivariate conditional copula cdf }
\item{tau2par}{function for maping Kendall's tau to the bivariate copula parameter}
}
\value{A list containing the following components:

\item{minimum}{the value of the estimated minimum of the negative log-likelihood}
\item{estimate}{the MLE}
\item{gradient}{the gradient at the estimated minimum of of the negative log-likelihood}
\item{hessian}{the hessian at the estimated minimum of the negative log-likelihood}
\item{code}{an integer indicating why the optimization process terminated}
\item{iterations}{the number of iterations performed}

For more details see  \code{\link{nlm}}
}
\references{
Nikoloulopoulos, A.K. (2016) Hybrid copula mixed models for combining case-control and cohort studies in meta-analysis  of diagnostic tests \emph{ArXiv e-prints}, arXiv:1604.05456. 
}
\seealso{
\code{\link{VineCopulaREMADA}}, \code{\link{CopulaREMADA}}
}
\examples{
# simulate the data from N=25 cohort studies
N=25
p=c(0.8,0.7,0.4)
g=c(0.1,0.1,0.05)
taus=c(-0.5,-0.3,-0.0001)
qcondcop12=qcondcop23=qcondcop13=qcondcln90
tau2par12=tau2par23=tau2par13=tau2par.cln90
simdat1=rVineCopulaREMADA.beta(N,p,g,taus,qcondcop12,qcondcop13,qcondcop23,
             tau2par12,tau2par13,tau2par23)
# simulate data from the N=25 case-control studies
tau=0.5
p=p[-3]
g=g[-3]
simdat2=rCopulaREMADA.beta(N,p,g,tau,rcln,tau2par.cln)
# combine the data
TP=c(simdat1$TP,simdat2$TP)
TN=c(simdat1$TN,simdat2$TN)
FP=c(simdat1$FP,simdat2$FP)
FN=c(simdat1$FN,simdat2$FN)
type=rep(c(1,2),each=N)

# fit the hybrid copula mixed model
nq=21
gl=gauss.quad.prob(nq,"uniform")
mgrid1<- meshgrid(gl$n,gl$n,gl$n,nargout=3)
mgrid2<- meshgrid(gl$n,gl$n)

perm=1
qcond=qcondcln
tau2par=tau2par.cln
# est=hybridCopulaREMADA.beta(TP,FN,FP,TN,type,perm,gl,mgrid1,mgrid2,
# qcondcop12,qcondcop13,tau2par12,tau2par13,qcond,tau2par)

}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

