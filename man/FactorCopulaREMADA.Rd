\name{FactorCopulaREMADA}
\Rdversion{1.1}
\alias{FactorCopulaREMADA}
\alias{FactorCopulaREMADA.norm}
\alias{FactorCopulaREMADA.beta}


\title{
Maximum likelihood estimation  of 1-factor copula mixed models for joint meta-analysis of \eqn{T} diagnostic tests
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
FactorCopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid,qcond1,tau2par1,qcond2,tau2par2)
                               
FactorCopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid,qcond1,tau2par1,qcond2,tau2par2)
}
\arguments{
\item{TP}{an \eqn{n\times T} matrix where \eqn{n} is the number of studies. Column \eqn{j} has  the number of true positives for test \eqn{j} for \eqn{j=1\ldots T}}
\item{FN}{an \eqn{n\times T} matrix where \eqn{n} is the number of studies. Column \eqn{j} has the number of false negatives Column \eqn{j} has  the number of true positives for test \eqn{j} for \eqn{j=1\ldots T}}
\item{FP}{an \eqn{n\times T} matrix where \eqn{n} is the number of studies. Column \eqn{j} has the number of false positives Column \eqn{j} has  the number of true positives for test \eqn{j} for \eqn{j=1\ldots T}}
\item{TN}{an \eqn{n\times T} matrix where \eqn{n} is the number of studies. Column \eqn{j} has the number of true negatives Column \eqn{j} has  the number of true positives for test \eqn{j} for \eqn{j=1\ldots T}}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{
a list containing two matrices with the rows of the output matrix x are copies of the vector \code{gl$nodes}; columns of the output matrix y are copies of the vector \code{gl$nodes}}
\item{qcond1}{function for the inverse conditional copula cdfs that link the factor with the latent sensitivities}
\item{tau2par1}{function for maping Kendall's tau  to copula parameter at the copulas that link the factor with the latent sensitivities}
\item{qcond2}{function for the inverse conditional copula cdfs that link the factor with the latent specificities}
\item{tau2par2}{function for maping Kendall's tau  to copula parameter at the copulas that link the factor with the latent specificities}
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
Nikoloulopoulos, A.K. (2022) An one-factor copula mixed model for joint meta-analysis of multiple diagnostic tests. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)},  185 (3), 1398--1423. \doi{10.1111/rssa.12838}.
}

\examples{
\donttest{
data(arthritis)
attach(arthritis)
TP=cbind(TP1,TP2)
TN=cbind(TN1,TN2)
FP=cbind(FP1,FP2)
FN=cbind(FN1,FN2)


nq=25
gl=gauss.quad.prob(nq,"uniform")
mgrid=meshgrid(gl$n,gl$n)
qcond1=qcondcln
qcond2=qcondcln270
tau2par1=tau2par.cln
tau2par2=tau2par.cln270

out=FactorCopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid,qcond1,tau2par1,qcond2,tau2par2)
se=sqrt(diag(solve(out$hessian)))

detach(arthritis)
}}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

