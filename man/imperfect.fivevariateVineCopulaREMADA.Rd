\name{imperfect.fivariateVineCopulaREMADA}
\Rdversion{1.1}
\alias{imperfect.fivevariateVineCopulaREMADA.norm.comprehensive}
\alias{imperfect.fivevariateVineCopulaREMADA.beta.comprehensive}



\title{
Maximum likelihood estimation of 5-variate 1-truncated D-vine copula mixed models for meta-analysis of diagnostic accuracy studies without a gold standard 
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
imperfect.fivevariateVineCopulaREMADA.norm.comprehensive(y11,y10,y01,y00,
                               gl,mgrid,qcond,tau2par,start)
                               
imperfect.fivevariateVineCopulaREMADA.beta.comprehensive(y11,y10,y01,y00,
                               gl,mgrid,qcond,tau2par,start)


}
\arguments{
\item{y11}{the number of the test results where the index test outcome is positive and the  reference test outcome is positive}
\item{y10}{the number of the test results where the index test outcome is positive and the  reference test outcome is negative}
\item{y01}{the number of the test results where the index test outcome is negative and the  reference test outcome is positive}
\item{y00}{the number of the test results where the index test outcome is negative and the  reference test outcome is negative}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{a list containing five-dimensional arrays.  Replicates of the quadrature points  that produce a 5-dimensional full grid}
\item{qcond}{function for the inverse of conditional copula cdf;  choices are \code{qconbvn} and \code{qcondfrk}}
\item{tau2par}{function for maping Kendall's tau to copula parameter; choices are \code{tau2par.bvn} and \code{tau2par.frk}}
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
\donttest{
data(Pap)
attach(Pap)

nq=15
gl=gauss.quad.prob(nq,"uniform")
data(mgrid5d15)
mgrid=mgrid5d15

tau2par=tau2par.bvn
qcond=qcondbvn

start=c(rep(0.6,5),rep(0.5,5),rep(0.01,4))
est.norm=imperfect.fivevariateVineCopulaREMADA.norm.comprehensive(y11,y10,y01,
y00,gl,mgrid,qcond,tau2par,start)

detach(Pap)
}}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

