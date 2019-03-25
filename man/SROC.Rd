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
SROC.norm(param,dcop,qcondcop,tau2par,TP,FN,FP,TN,
          points=TRUE,curves=TRUE,
          NEP=rep(0,length(TP)),NEN=rep(0,length(TP)))
SROC.beta(param,dcop,qcondcop,tau2par,TP,FN,FP,TN,
          points=TRUE,curves=TRUE,
          NEP=rep(0,length(TP)),NEN=rep(0,length(TP)))
SROC(param.beta,param.normal,TP,FN,FP,TN,
          NEP=rep(0,length(TP)),NEN=rep(0,length(TP)))
}
\arguments{
\item{param}{A vector with the sensitivities, specifities, variabilities and Kendall's tau value (the latter only for \code{\link{SROC.norm}} and \code{\link{SROC.beta}})}
\item{param.beta}{A vector with the sensitivity, specifity and variabilities of the countermonotonic CopulaREMADA with beta margins}
\item{param.normal}{A vector with the sensitivity, specifity and variabilities of the countermonotonic CopulaREMADA with normal margins}
\item{dcop}{function for copula density}
\item{qcondcop}{function for the inverse of conditional copula cdf}
\item{tau2par}{function for maping Kendall's tau to copula parameter}
\item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
\item{points}{logical: print individual studies}
\item{curves}{logical: print quantile regression curves}
\item{NEP}{the number of non-evaluable positives in the presence of non-evaluable subjects}
\item{NEN}{the number of non-evaluable negatives in the presence of non-evaluable subjects}
}
\value{Summary receiver operating characteristic   curves
}
\references{
Nikoloulopoulos, A.K.  (2015) A mixed effect model for bivariate meta-analysis of diagnostic test accuracy studies using a copula representation of the random effects distribution. \emph{Statistics in Medicine}, \bold{34}, 3842--3865. \doi{10.1002/sim.6595}.
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

data(MK2016)
attach(MK2016)
p=c(0.898745016,0.766105342,0.059168715,0.109217888)
g=c(0.090270947,0.079469009,0.367463579,0.154976269)
taus=c(0.82050793,-0.51867629,0.26457961)
SROC.beta(c(p[1:2],g[1:2],taus[1]),
          dcln180,qcondcln180,tau2par.cln180,
          TP,FN,FP,TN,points=TRUE,curves=TRUE,NEP,NEN)
detach(MK2016)


}
\keyword{graphics}
