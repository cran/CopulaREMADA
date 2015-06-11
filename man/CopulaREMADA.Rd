\name{CopulaREMADA}
\Rdversion{1.1}
\alias{CopulaREMADA}
\alias{CopulaREMADA.norm}
\alias{CopulaREMADA.beta}
\alias{countermonotonicCopulaREMADA}
\alias{countermonotonicCopulaREMADA.norm}
\alias{countermonotonicCopulaREMADA.beta}
\title{
Maximum likelhood estimation  for copula mixed models for diagnostic test accurracy studies
}
\description{
For  copula mixed models for diagnostic test accuracy studies numerical evaluation of the MLE is easily done with the following steps:


1. Calculate Gauss-Legendre  quadrature points 
{gl$nodes} and weights {gl$weights}.

2. Convert from independent uniform quadrature points  to dependent uniform quadrature points  that have distribution 'cop'.
The inverse of the conditional distribution \code{qcondcop} corresponding to the copula 'cop' is used  to achieve this.

3. Numerically evaluate the joint probability mass function with the bivariate integral
in a double sum. 


With Gauss-Legendre quadrature, the same nodes and weights
are used for different functions;
this helps in yielding smooth numerical derivatives for numerical optimization via quasi-Newton.
Our comparisons show that \eqn{n_q=15} is adequate with good precision to at least at four decimal places. 
}
\usage{
CopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid,qcond,tau2par)
CopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid,qcond,tau2par) 
countermonotonicCopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid) 
countermonotonicCopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid) 
}
\arguments{
\item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{
a list containing two matrices with the rows of the output matrix x are copies of the vector \code{gl$nodes}; columns of the output matrix y are copies of the vector \code{gl$nodes}. For more details see also \code{\link{meshgrid}}}
\item{qcond}{function for conditional copula cdf}
\item{tau2par}{function for maping Kendall's tau to copula parameter}
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
Nikoloulopoulos, A.K. (2015a) A mixed effect model for bivariate meta-analysis of diagnostic test accuracy studies using  a copula representation of the random effects distribution, \emph{ArXiv e-prints}. 
}
\seealso{
\code{\link{rCopulaREMADA}}
}
\examples{
nq=15
gl=gauss.quad.prob(nq,"uniform")
mgrid<- meshgrid(gl$n,gl$n)

data(LAG)
attach(LAG)
c270est.b=CopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid,qcondcln270,tau2par.cln270)
detach(LAG)

data(MRI)
attach(MRI)
c270est.n=CopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid,qcondcln270,tau2par.cln270)
detach(MRI)

data(CT)
attach(CT)
est.n=countermonotonicCopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid)
est.b=countermonotonicCopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid)
detach(CT)
}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

