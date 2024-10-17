\name{imperfectREMADA}
\Rdversion{1.1}
\alias{imperfectREMADA}
\alias{imperfectREMADA.norm}
\alias{imperfectREMADA.beta}


\title{
Maximum likelihood estimation of univariate mixed models for meta-analysis of diagnostic accuracy studies without a gold standard 
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
imperfectREMADA.norm(y11,y10,y01,y00,gl,select.random,start)
                               
imperfectREMADA.beta(y11,y10,y01,y00,gl,select.random,start)
}
\arguments{
\item{y11}{the number of the test results where the index test outcome is positive and the  reference test outcome is positive}
\item{y10}{the number of the test results where the index test outcome is positive and the  reference test outcome is negative}
\item{y01}{the number of the test results where the index test outcome is negative and the  reference test outcome is positive}
\item{y00}{the number of the test results where the index test outcome is negative and the  reference test outcome is negative}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{select.random}{a scalar from 1 to 5 that indicates the random effect}
\item{start}{starting values for the parameters}
}
\value{A list containing the following components:
\item{LogLikelihood}{the maximized log-likelihood}
\item{Estimates}{the MLE}
\item{SE}{the standard errors}
}
\references{
Nikoloulopoulos, A.K. (2024) Vine copula mixed models for meta-analysis of diagnostic accuracy studies without a gold standard. Submitted.
}

\examples{
data(Pap)
attach(Pap)

nq=30
gl=gauss.quad.prob(nq,"uniform")
start=c(rep(0.6,5),0.5)
select.random=1 
est.norm=imperfectREMADA.norm(y11,y10,y01,y00,gl,
select.random,start)
start=c(rep(0.6,5),0.1)
est.beta=imperfectREMADA.beta(y11,y10,y01,y00,gl,
select.random,start)

detach(Pap)
}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

