\name{quadVineCopulaREMADA}
\Rdversion{1.1}
\alias{quadVineCopulaREMADA}
\alias{quadVineCopulaREMADA.norm}
\alias{quadVineCopulaREMADA.beta}
\alias{quadVineCopulaREMADA.norm.beta}


\title{
Maximum likelihood estimation  of quadrivariate  D-vine copula mixed models for joint meta-analysis and comparison of  two diagnostic tests
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
quadVineCopulaREMADA.norm(TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
                               gl,mgrid,qcond1,pcond1,tau2par1,
                               qcond2,pcond2,tau2par2)
                               
quadVineCopulaREMADA.beta(TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
                               gl,mgrid,qcond1,pcond1,tau2par1,
                               qcond2,pcond2,tau2par2) 
quadVineCopulaREMADA.norm.beta(TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
                               gl,mgrid,qcond1,pcond1,tau2par1,
                               qcond2,pcond2,tau2par2)              
}
\arguments{
\item{TP1}{the number of true positives for test 1}
\item{FN1}{the number of false negatives for test 1}
\item{FP1}{the number of false positives for test 1}
\item{TN1}{the number of true negatives for test 1}
\item{TP2}{the number of true positives for test 2}
\item{FN2}{the number of false negatives for test 2}
\item{FP2}{the number of false positives for test 2}
\item{TN2}{the number of true negatives for test 2}

\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{a list containing four-dimensional arrays.  Replicates of the quadrature points  that produce a 4-dimensional full grid}
\item{qcond1}{function for the inverse conditional copula cdf at the (1,2) bivariate margin}
\item{pcond1}{function for the conditional copula cdf at the (1,2) bivariate margin}
\item{tau2par1}{function for maping Kendall's tau at the (1,2) bivariate margin to copula parameter}
\item{qcond2}{function for the inverse conditional copula cdf at the (3,4) bivariate margin}
\item{pcond2}{function for the conditional copula cdf at the (3,4) bivariate margin}
\item{tau2par2}{function for maping Kendall's tau at the (3,4) bivariate margin to copula parameter}
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
Nikoloulopoulos, A.K. (2019) A D-vine copula mixed model for joint meta-analysis and comparison of  diagnostic tests.  \emph{Statistical Methods in Medical Research}, 28(10-11):3286--3300. \doi{10.1177/0962280218796685}. 
}

\examples{
\donttest{
nq=15
gl=gauss.quad.prob(nq,"uniform")
data(mgrid15)

data(arthritis)
attach(arthritis)

qcond1=qcondcln270
pcond1=pcondcln270
tau2par1=tau2par.cln270

qcond2=qcondfrk
pcond2=pcondfrk
tau2par2=tau2par.frk

out<-quadVineCopulaREMADA.norm(TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
gl,mgrid15,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)

detach(arthritis)
}}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

