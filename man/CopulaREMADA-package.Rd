\name{CopulaREMADA-package}
\alias{CopulaREMADA-package}
\docType{package}
\title{Copula Mixed  Models for Multivariate Meta-Analysis of Diagnostic Test Accuracy Studies}
\description{
 Fits copula mixed models for multivariate meta-analysis of diagnostic test accuracy studies proposed in Nikoloulopoulos (2015, 2017, 2018a, 2018b, 2019, 2020a, 2020b, 2020c). 
}
\details{
This package contains \R functions to implement: 
\itemize{
\item The copula mixed model for meta-analysis of diagnostic test accuracy studies and produce SROC curves and summary operating points (a pair of average sensitivity and specificity) with a confidence region and a predictive region (Nikoloulopoulos, 2015, 2018a). All the analyses presented in Section 7 of Nikoloulopoulos (2015) are given as code examples in the package;


\item The vine copula mixed model for meta-analysis of diagnostic test accuracy studies accounting for disease prevalence and non-evaluable subjects (Nikoloulopoulos, 2017, 2020a); 


\item The hybrid vine copula mixed model for meta-analysis of diagnostic test accuracy case-control and cohort studies (Nikoloulopoulos, 2018b);


\item The D-vine copula mixed model for meta-analysis and comparison of two diagnostic tests (Nikoloulopoulos, 2019). 

\item The multinomial quadrivariate D-vine copula mixed model for diagnostic studies meta-analysis accounting for non-evaluable subjects (Nikoloulopoulos, 2020b). 

\item The one-factor copula mixed model for joint meta-analysis of multiple diagnostic tests (Nikoloulopoulos, 2020c). 
}



}
\author{
Aristidis K. Nikoloulopoulos.
}
\references{

Nikoloulopoulos, A.K.  (2015) A mixed effect model for bivariate meta-analysis of diagnostic test accuracy studies using a copula representation of the random effects distribution. \emph{Statistics in Medicine}, \bold{34}, 3842--3865. \doi{10.1002/sim.6595}.

Nikoloulopoulos, A.K.  (2017) A vine copula mixed effect model for trivariate meta-analysis of diagnostic test accuracy studies accounting for disease prevalence. \emph{Statistical Methods in Medical Research}, \bold{26}, 2270--2286. \doi{10.1177/0962280215596769}.

Nikoloulopoulos, A.K. (2018a) On composite likelihood in bivariate meta-analysis of diagnostic test accuracy studies.  \emph{AStA Advances in Statistical Analysis}, \bold{102}, 211--227. \doi{10.1007/s10182-017-0299-y}. 

Nikoloulopoulos, A.K. (2018b) Hybrid copula mixed models for combining case-control and cohort studies in meta-analysis of diagnostic tests. \emph{Statistical Methods in Medical Research}, \bold{27}, 2540--2553. \doi{10.1177/0962280216682376}. 

Nikoloulopoulos, A.K.  (2019) A D-vine copula mixed model for joint meta-analysis and comparison of  diagnostic tests.  \emph{Statistical Methods in Medical Research},  28(10-11):3286--3300. \doi{10.1177/0962280218796685}. 

Nikoloulopoulos, A.K.  (2020a) An extended trivariate vine copula mixed model for meta-analysis of diagnostic studies in the presence of non-evaluable outcomes.  \emph{The International Journal of Biostatistics}, 16(2). \doi{10.1515/ijb-2019-0107}.


Nikoloulopoulos, A.K. (2020b) A multinomial quadrivariate D-vine copula mixed model for diagnostic studies meta-analysis in the presence of non-evaluable subjects.
\emph{Statistical Methods in Medical Research},  29 (10), 2988--3005. \doi{10.1177/0962280220913898}. 

Nikoloulopoulos, A.K. (2020c) An one-factor copula mixed model for joint meta-analysis of multiple diagnostic tests. 
\emph{ArXiv e-prints}, arXiv:2006.09278. \url{https://arxiv.org/abs/2006.09278}.


}
\keyword{package}
