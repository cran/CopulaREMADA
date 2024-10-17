\name{vine.vuong}
\Rdversion{1.1}
\alias{vine.vuong}
\alias{vine.vuong.norm}
\alias{vine.vuong.beta}
\alias{tvine.vuong.norm}
\alias{tvine.vuong.beta}
\alias{tvine2.vuong.norm}
\alias{tvine2.vuong.beta}

\title{
Vuong's test for the comparison of non-nested vine copula mixed models for diagnostic test accuaracy studies
}
\description{
Vuong (1989)'s test for the comparison of non-nested vine copula mixed models for diagnostic test accuaracy studies.  
It shows if a vine copula mixed model provides better fit than the standard GLMM. We compute the Vuong's test with Model 1 being the vine copula mixed model with BVN copula and normal margins, i.e., the standard GLMM.
}
\usage{
vine.vuong.beta(qcondcop12,qcondcop13,qcondcop23,
tau2par12,tau2par13,tau2par23,param1,param2,TP,FN,FP,TN,gl,mgrid,NEP,NEN)
vine.vuong.norm(qcondcop12,qcondcop13,qcondcop23,
tau2par12,tau2par13,tau2par23,param1,param2,TP,FN,FP,TN,gl,mgrid,NEP,NEN) 
tvine.vuong.beta(qcondcop12,qcondcop13,
tau2par12,tau2par13,param1,param2,TP,FN,FP,TN,gl,mgrid,NEP,NEN)
tvine.vuong.norm(qcondcop12,qcondcop13,
tau2par12,tau2par13,param1,param2,TP,FN,FP,TN,gl,mgrid,NEP,NEN) 
tvine2.vuong.beta(qcondcop12,qcondcop13,
tau2par12,tau2par13,param1,param2,TP,FN,FP,TN,gl,mgrid,NEP,NEN)
tvine2.vuong.norm(qcondcop12,qcondcop13,
tau2par12,tau2par13,param1,param2,TP,FN,FP,TN,gl,mgrid,NEP,NEN) 
}
\arguments{
\item{qcondcop12}{function for the inverse of conditional copula cdf at the (1,2) bivariate margin for Model 2}
\item{qcondcop13}{function for the inverse of conditional copula cdf at the (1,3) bivariate margin for Model 2}
\item{qcondcop23}{function for the inverse of conditional copula cdf at the (2,3|1) bivariate margin for Model 2}
\item{tau2par12}{function for maping Kendall's tau at the (1,2) bivariate margin to copula parameter for Model 2}
\item{tau2par13}{function for maping Kendall's tau at the (1,3) bivariate margin to copula parameter for Model 2}
\item{tau2par23}{function for maping Kendall's tau at the (2,3|1) bivariate margin to the conditional copula parameter for Model 2}
\item{param1}{parameters for the Model 1. i.e., the GLMM}
\item{param2}{parameters for the Model 2}
\item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{
a list containing three-dimensional arrays}
\item{NEP}{the number of non-evaluable positives in the presence of non-evaluable subjects}
\item{NEN}{the number of non-evaluable negatives in the presence of non-evaluable subjects}
}
\value{A list containing the following components:

\item{z}{the test statistic}
\item{p-value}{the \eqn{p}-value}
}
\references{
Nikoloulopoulos, A.K.  (2017) A vine copula mixed effect model for trivariate meta-analysis of diagnostic test accuracy studies accounting for disease prevalence. \emph{Statistical Methods in Medical Research}, \bold{26}, 2270--2286. \doi{10.1177/0962280215596769}.

Nikoloulopoulos, A.K. (2020) An extended trivariate vine copula mixed model for meta-analysis of diagnostic studies in the presence of non-evaluable outcomes.  \emph{The International Journal of Biostatistics}, 16(2). \doi{10.1515/ijb-2019-0107}.

Vuong Q.H. (1989) Likelihood ratio tests for model selection and non-nested hypotheses. \emph{Econometrica}, \bold{57}, 307--333.
}
\seealso{
\code{\link{CopulaREMADA}}
}
\examples{
nq=15
gl=gauss.quad.prob(nq,"uniform")
mgrid=meshgrid(gl$n,gl$n,gl$n,nargout=3)

data(betaDG)
attach(betaDG)
#nest.n2=VineCopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid,
#qcondbvn,qcondbvn,qcondbvn,
#tau2par.bvn,tau2par.bvn,tau2par.bvn)
nest.n2.est= #nest.n2$e
c(0.87186926,  0.13696066,  0.70614956,  0.69152133,  
0.51780203,  0.70883558, -0.41354870,0.07701287, -0.12111253)
#c090est.b2=VineCopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid,
#qcondcln90,qcondcln,qcondcln90,tau2par.cln90,tau2par.cln,tau2par.cln90)
c090est.b2.est= #c090est.b2$e
c(0.85528463,  0.14667571,  0.68321231,  0.04897466,
0.02776290,  0.08561436, -0.34639172, 0.04621924, -0.21627977)
c090vuong.b2=vine.vuong.beta(qcondcln90,qcondcln,qcondcln90,
tau2par.cln90,tau2par.cln,tau2par.cln90,
nest.n2.est,c090est.b2.est,TP,FN,FP,TN,gl,mgrid)
c090vuong.b2
detach(betaDG)
}
\keyword{maximum likelihood}

