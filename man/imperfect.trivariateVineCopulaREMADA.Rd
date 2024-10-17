\name{imperfect.trivariateVineCopulaREMADA}
\Rdversion{1.1}
\alias{imperfect.trivariateVineCopulaREMADA}
\alias{imperfect.trivariateVineCopulaREMADA.norm}
\alias{imperfect.trivariateVineCopulaREMADA.beta}
\alias{imperfect.trivariateVineCopulaREMADA.norm.comprehensive}
\alias{imperfect.trivariateVineCopulaREMADA.beta.comprehensive}



\title{
Maximum likelihood estimation of trivariate 1-truncated D-vine copula mixed models for meta-analysis of diagnostic accuracy studies without a gold standard 
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
imperfect.trivariateVineCopulaREMADA.norm.comprehensive(y11,
y10,y01,y00,gl,mgrid,qcond,tau2par,select.random,start)
                               
imperfect.trivariateVineCopulaREMADA.norm(y11,
y10,y01,y00,gl,mgrid,qcond1,tau2par1,qcond2,tau2par2,select.random,start)

imperfect.trivariateVineCopulaREMADA.beta.comprehensive(y11,
y10,y01,y00,gl,mgrid,qcond,tau2par,select.random,start)
                               
imperfect.trivariateVineCopulaREMADA.beta(y11,
y10,y01,y00,gl,mgrid,qcond1,tau2par1,qcond2,tau2par2,select.random,start)

}
\arguments{
\item{y11}{the number of the test results where the index test outcome is positive and the  reference test outcome is positive}
\item{y10}{the number of the test results where the index test outcome is positive and the  reference test outcome is negative}
\item{y01}{the number of the test results where the index test outcome is negative and the  reference test outcome is positive}
\item{y00}{the number of the test results where the index test outcome is negative and the  reference test outcome is negative}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{a list containing three-dimensional arrays.  Replicates of the quadrature points  that produce a 3-dimensional full grid}
\item{select.random}{vector \eqn{(t_{1},t_{2},t_3)}, where \eqn{1\leq t_1<t_2<t_3\leq 5}, that indicates the random effects}
\item{qcond}{function for the inverse of conditional copula cdf;  choices are \code{qconbvn} and \code{qcondfrk}}
\item{tau2par}{function for maping Kendall's tau to copula parameter; choices are \code{tau2par.bvn} and \code{tau2par.frk}}
\item{qcond1}{function for the inverse of conditional copula cdf for the \eqn{(t_{1},t_{2})} bivariate margin; choices are \code{qcondcln} and \code{qcondcln180}}
\item{tau2par1}{function for maping Kendall's tau to copula parameter for the \eqn{(t_{1},t_{2})} bivariate margin; choices are \code{tau2par.cln} and \code{tau2par.cln180}}
\item{qcond2}{function for the inverse of conditional copula cdf for the \eqn{(t_{2},t_{3})} bivariate margin; choices are \code{qcondcln90} and \code{qcondcln270}}
\item{tau2par2}{function for maping Kendall's tau to copula parameter for the \eqn{(t_{2},t_{3})} bivariate margin; choices are \code{tau2par.cln90} and \code{tau2par.cln270}}
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

nq=30
gl=gauss.quad.prob(nq,"uniform")
mgrid<- meshgrid(gl$n,gl$n,gl$n,nargout=3)

tau2par=tau2par.bvn
qcond=qcondbvn

select.random=c(1,2,4)
start=c(rep(0.6,5),rep(0.5,3),c(0.01,-0.01))
est.norm=imperfect.trivariateVineCopulaREMADA.norm.comprehensive(y11,y10,
y01,y00,gl,mgrid,qcond,tau2par,select.random,start)

tau2par1=tau2par.cln180
qcond1=qcondcln180
tau2par2=tau2par.cln270
qcond2=qcondcln270
est.norm.cln=imperfect.trivariateVineCopulaREMADA.norm(y11,y10,y01,
y00,gl,mgrid,qcond1,tau2par1,qcond2,tau2par2,select.random,start)

start=c(rep(0.6,5),rep(0.05,3),c(0.1,-0.1))
est.beta=imperfect.trivariateVineCopulaREMADA.beta.comprehensive(y11,y10,y01,y00,
gl,mgrid,qcond,tau2par,select.random,start)

est.beta.cln=imperfect.trivariateVineCopulaREMADA.beta(y11,y10,y01,y00,
gl,mgrid,qcond1,tau2par1,qcond2,tau2par2,select.random,start)

detach(Pap)
}}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

