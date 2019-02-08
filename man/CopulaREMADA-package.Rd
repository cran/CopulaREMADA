\name{CopulaREMADA-package}
\alias{CopulaREMADA-package}
\docType{package}
\title{Copula Mixed Effect Models for Multivariate Meta-Analysis of Diagnostic Test Accuracy Studies}
\description{
 Fits copula mixed effect models for multivariate meta-analysis of diagnostic test accuracy studies proposed in Nikoloulopoulos (2015, 2017, 2018a, 2018b, 2018c). 
}
\details{
This package contains \R functions to implement: 
\itemize{
\item The copula mixed model for meta-analysis of diagnostic test accuracy studies and produce SROC curves and summary operating points (a pair of average sensitivity and specificity) with a confidence region and a predictive region (Nikoloulopoulos, 2015, 2018a). All the analyses presented in Section 7 of Nikoloulopoulos (2015) are given as code examples in the package;


\item The vine copula mixed model for meta-analysis of diagnostic test accuracy studies accounting for disease prevalence and non-evaluable subjects (Nikoloulopoulos, 2017, 2018b); 


\item The hybrid vine copula mixed model for meta-analysis of diagnostic test accuracy case-control and cohort studies (Nikoloulopoulos, 2018c);


\item The D-vine copula mixed model for meta-analysis and comparison of two diagnostic tests (Nikoloulopoulos, 2018d).  
}



}
\author{
Aristidis K. Nikoloulopoulos.
}
\references{

Nikoloulopoulos, A.K.  (2015) A mixed effect model for bivariate meta-analysis of diagnostic test accuracy studies using a copula representation of the random effects distribution. \emph{Statistics in Medicine}, \bold{34}, 3842--3865. \doi{10.1002/sim.6595}.

Nikoloulopoulos, A.K.  (2017) A vine copula mixed effect model for trivariate meta-analysis of diagnostic test accuracy studies accounting for disease prevalence. \emph{Statistical Methods in Medical Research}, \bold{26}, 2270--2286. \doi{10.1177/0962280215596769}.

Nikoloulopoulos, A.K. (2018a) On composite likelihood in bivariate meta-analysis of diagnostic test accuracy studies.  \emph{AStA Advances in Statistical Analysis}, \bold{102}, 211--227. \doi{10.1007/s10182-017-0299-y}. 

Nikoloulopoulos, A.K.  (2018b) A vine copula mixed  model for trivariate meta-analysis of diagnostic studies accounting for disease prevalence and non-evaluable subjects.  \emph{ArXiv e-prints}, arXiv:1812.03685. \url{https://arxiv.org/abs/1812.03685}.

Nikoloulopoulos, A.K. (2018c) Hybrid copula mixed models for combining case-control and cohort studies in meta-analysis of diagnostic tests. \emph{Statistical Methods in Medical Research}, \bold{27}, 2540--2553. \doi{10.1177/0962280216682376}. 

Nikoloulopoulos, A.K.  (2018d) A D-vine copula mixed model for joint meta-analysis and comparison of  diagnostic tests.  \emph{Statistical Methods in Medical Research}, in press. \doi{10.1177/0962280218796685}. 




}
\keyword{package}
