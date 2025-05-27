\name{imperfectCopulaREMADA}
\Rdversion{1.1}
\alias{imperfectCopulaREMADA}
\alias{imperfectCopulaREMADA.norm}
\alias{imperfectCopulaREMADA.beta}


\title{
Maximum likelihood estimation of bivariate copula mixed models for meta-analysis of diagnostic accuracy studies without a gold standard 
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
imperfectCopulaREMADA.norm(y11,y10,y01,y00,
gl,mgrid,qcond,tau2par,select.random,start)
                               
imperfectCopulaREMADA.beta(y11,y10,y01,y00,
gl,mgrid,qcond,tau2par,select.random,start)

}
\arguments{
\item{y11}{the number of the test results where the index test outcome is positive and the  reference test outcome is positive}
\item{y10}{the number of the test results where the index test outcome is positive and the  reference test outcome is negative}
\item{y01}{the number of the test results where the index test outcome is negative and the  reference test outcome is positive}
\item{y00}{the number of the test results where the index test outcome is negative and the  reference test outcome is negative}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{
a list containing two matrices with the rows of the output matrix x are copies of the vector \code{gl$nodes}; columns of the output matrix y are copies of the vector \code{gl$nodes}}
\item{select.random}{vector \eqn{(t_{1},t_{2})}, where \eqn{1\leq t_1<t_2\leq 5}, that indicates the random effects}
\item{qcond}{function for the inverse of conditional copula cdf}
\item{tau2par}{function for maping Kendall's tau to copula parameter}
\item{start}{starting values for the parameters}
}
\value{A list containing the following components:

\item{LogLikelihood}{the maximized log-likelihood}
\item{Estimates}{the MLE}
\item{SE}{the standard errors}
}
\references{
Nikoloulopoulos, A.K. (2025) Vine copula mixed models for meta-analysis of diagnostic accuracy studies without a gold standard. \emph{Biometrics}, \bold{81}(2), ujaf037. \doi{10.1093/biomtc/ujaf037}.
}

\examples{
data(Pap)
attach(Pap)

nq=30
gl=gauss.quad.prob(nq,"uniform")
mgrid<- meshgrid(gl$n,gl$n)

tau2par=tau2par.bvn
qcond=qcondbvn

select.random=c(1,2)

start=c(rep(0.6,5),rep(0.5,2),-0.1)
est.norm=imperfectCopulaREMADA.norm(y11,y10,y01,y00,gl,mgrid,
qcond,tau2par,select.random,start)

detach(Pap)
}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

