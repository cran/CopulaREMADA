\name{vuong}
\Rdversion{1.1}
\alias{vuong}
\alias{countermonotonicity.vuong}
\alias{vuong.norm}
\alias{vuong.beta}

\title{
Vuong's test for the comparison of non-nested copula mixed models for diagnostic test accuaracy studies
}
\description{
Vuong (1989)'s test for the comparison of non-nested copula mixed models for diagnostic test accuaracy studies.  
It shows if a copula mixed model provides better fit than the standard GLMM. We compute the Vuong's test with Model 1 being the copula mixed model with BVN copula and normal margins, i.e., the standard GLMM.
}
\usage{
vuong.norm(qcond,tau2par,param1,param2,TP,FN,FP,TN,gl,mgrid)
vuong.beta(qcond,tau2par,param1,param2,TP,FN,FP,TN,gl,mgrid)
countermonotonicity.vuong(param1,param2,TP,FN,FP,TN,gl,mgrid)  
}
\arguments{
\item{qcond}{function for conditional copula cdf for Model 2}
\item{tau2par}{function for maping Kendall's tau to copula parameter for Model 2}
\item{param1}{parameters for the Model 1. i.e., the GLMM}
\item{param2}{parameters for the Model 2}
\item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{
a list containing two matrices with the rows of the output matrix X are copies of the vector \code{gl$nodes}; columns of the output matrix Y are copies of the vector \code{gl$nodes}. For more details see \code{\link{meshgrid}}}
}
\value{A list containing the following components:

\item{nu}{the test statistic}
\item{p-value}{the \eqn{p}-value}
}
\references{
Nikoloulopoulos, A.K. (2015) A mixed effect model for bivariate meta-analysis of diagnostic test accuracy studies using  a copula representation of the random effects distribution, \emph{ArXiv e-prints}. 

Vuong Q.H. (1989) Likelihood ratio tests for model selection and non-nested hypotheses. \emph{Econometrica} 57(2), 307--333.
}
\seealso{
\code{\link{CopulaREMADA}}
}
\examples{
library(statmod)
nq=15
gl=gauss.quad.prob(nq,"uniform")
mgrid<- meshgrid(gl$n,gl$n)

data(MRI)
attach(MRI)
c270est.b=CopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid,qcondcln270,tau2par.cln270)
nest.n=CopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid,qcondbvn,tau2par.bvn)
c90est.n=CopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid,qcondcln90,tau2par.cln90)
vuong.beta(qcondcln270,tau2par.cln270,nest.n$e,c270est.b$e,TP,FN,FP,TN,gl,mgrid)
vuong.norm(qcondcln90,tau2par.cln90,nest.n$e,c90est.n$e,TP,FN,FP,TN,gl,mgrid)
detach(MRI)

data(CT)
attach(CT)
est.n=countermonotonicCopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid)
est.b=countermonotonicCopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid)
countermonotonicity.vuong(est.n$e,est.b$e,TP,FN,FP,TN,gl,mgrid)
detach(CT)
}
\keyword{maximum likelihood}

