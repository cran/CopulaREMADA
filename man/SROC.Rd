\name{SROC}
\Rdversion{1.1}
\alias{SROC.norm}
\alias{SROC.beta}
\alias{SROC}
\title{
Summary receiver operating characteristic   curves for copula mixed effect models for bivariate meta-analysis of diagnostic test accuracy studies
}
\description{
Summary receiver operating characteristic (SROC) curves are demonstrated for the proposed models through quantile regression techniques and different characterizations of the estimated bivariate random effects distribution  
}
\usage{
SROC.norm(param,dcop,qcondcop,tau2par,TP,FN,FP,TN,points=TRUE,curves=TRUE)
SROC.beta(param,dcop,qcondcop,tau2par,TP,FN,FP,TN,points=TRUE,curves=TRUE)
SROC(param.beta,param.normal,TP,FN,FP,TN)
}
\arguments{
\item{param}{A vector with the sensitivities, specifities, variabilities and Kendall's tau value (the latter only for \code{\link{SROC.norm}} and \code{\link{SROC.beta}})}
\item{param.beta}{A vector with the sensitivity, specifity and variabilities of the countermonotonic CopulaREMADA with beta margins}
\item{param.normal}{A vector with the sensitivity, specifity and variabilities of the countermonotonic CopulaREMADA with normal margins}
\item{dcop}{function for copula density}
\item{qcondcop}{function for conditional copula cdf}
\item{tau2par}{function for maping Kendall's tau to copula parameter}
\item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
\item{points}{logical: print individual studies}
\item{curves}{logical: print quantile regression curves}
}
\value{Summary receiver operating characteristic   curves
}
\references{
Nikoloulopoulos, A.K. (2015a) A mixed effect model for bivariate meta-analysis of diagnostic test accuracy studies using  a copula representation of the random effects distribution, \emph{ArXiv e-prints}. 
}

\seealso{
\code{\link{CopulaREMADA}}
\code{\link{rCopulaREMADA}}
}
\examples{
nq=15
gl=gauss.quad.prob(nq,"uniform")
mgrid<- meshgrid(gl$n,gl$n)

data(telomerase) 
attach(telomerase)
est.n=countermonotonicCopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid)
est.b=countermonotonicCopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid)
SROC(est.b$e,est.n$e,TP,FN,FP,TN)
detach(telomerase)

data(LAG)
attach(LAG)
c180est.b=CopulaREMADA.beta(TP,FN,FP,TN,gl,mgrid,qcondcln180,tau2par.cln180)
SROC.beta(c180est.b$e,dcln180,qcondcln180,tau2par.cln180,TP,FN,FP,TN)
detach(LAG)

data(MRI)
attach(MRI)
c270est.n=CopulaREMADA.norm(TP,FN,FP,TN,gl,mgrid,qcondcln270,tau2par.cln270)
SROC.norm(c270est.n$e,dcln270,qcondcln270,tau2par.cln270,TP,FN,FP,TN)
detach(MRI)

}
\keyword{graphics}
