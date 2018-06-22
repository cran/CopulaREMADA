\name{VineCopulaREMADA}
\Rdversion{1.1}
\alias{VineCopulaREMADA}
\alias{VineCopulaREMADA.norm}
\alias{VineCopulaREMADA.beta}
\alias{tVineCopulaREMADA}
\alias{tVineCopulaREMADA.norm}
\alias{tVineCopulaREMADA.beta}

\title{
Maximum likelhood estimation  for (truncated) vine copula mixed models for diagnostic test accurracy studies accounting for disease prevalence
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
VineCopulaREMADA.norm(TP,FN,FP,TN,perm,gl,mgrid,
                               qcondcop12,qcondcop13,qcondcop23,
                               tau2par12,tau2par13,tau2par23)
VineCopulaREMADA.beta(TP,FN,FP,TN,perm,gl,mgrid,
                               qcondcop12,qcondcop13,qcondcop23,
                               tau2par12,tau2par13,tau2par23) 
tVineCopulaREMADA.norm(TP,FN,FP,TN,perm,gl,mgrid,
                               qcondcop12,qcondcop13,
                               tau2par12,tau2par13)
tVineCopulaREMADA.beta(TP,FN,FP,TN,perm,gl,mgrid,
                               qcondcop12,qcondcop13,
                               tau2par12,tau2par13) 
}
\arguments{
\item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
\item{perm}{a scalar indicating the selected permutation of indices: 

1: Pilot variable is the number of TP. The bivariate margins are 12, 13, 23|1; 

2: Pilot variable is the number of TN. The bivariate margins are 23, 12, 13|2; 

3: Pilot variable is the TP+FN. The bivariate margins are 13, 23, 12|3; 

1:TP, 2:TN, 3:TP+FN}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{
a list containing three-dimensional arrays. For more details see \code{\link{meshgrid}}}
\item{qcondcop12}{function for the inverse of conditional copula cdf at the (1,2) bivariate margin}
\item{qcondcop13}{function for the inverse of conditional copula cdf at the (1,3) bivariate margin}
\item{qcondcop23}{function for the inverse of conditional copula cdf at the (2,3|1) bivariate margin}
\item{tau2par12}{function for maping Kendall's tau at the (1,2) bivariate margin to copula parameter}
\item{tau2par13}{function for maping Kendall's tau at the (1,3) bivariate margin to copula parameter}
\item{tau2par23}{function for maping Kendall's tau at the (2,3|1) bivariate margin to the conditional copula parameter}
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
Nikoloulopoulos, A.K.  (2017) A vine copula mixed effect model for trivariate meta-analysis of diagnostic test accuracy studies accounting for disease prevalence. \emph{Statistical Methods in Medical Research}, \bold{26}, 2270--2286. \doi{10.1177/0962280215596769}. 
}
\seealso{
\code{\link{rVineCopulaREMADA}}
}
\examples{
nq=15
gl=gauss.quad.prob(nq,"uniform")
mgrid=meshgrid(gl$n,gl$n,gl$n,nargout=3)

data(OGT)
attach(OGT)
#out=tVineCopulaREMADA.norm(TP,FN,FP,TN,1,gl,mgrid,
#qcondbvn,qcondbvn,tau2par.bvn,tau2par.bvn)
detach(OGT)
}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

