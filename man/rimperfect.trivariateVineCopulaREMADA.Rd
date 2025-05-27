\name{rimperfect.trivariateVineCopulaREMADA}
\Rdversion{1.1}
\alias{rimperfect.trivariateVineCopulaREMADA}
\alias{rimperfect.trivariateVineCopulaREMADA.norm}
\alias{rimperfect.trivariateVineCopulaREMADA.beta}
\title{
Simulation from  trivariate 1-truncated D-vine copula mixed models for meta-analysis of diagnostic accuracy studies without a gold standard 
}

\description{
Simulation from trivariate 1-truncated D-vine copula mixed models for meta-analysis of diagnostic accuracy studies without a gold standard
}

\usage{
rimperfect.trivariateVineCopulaREMADA.norm(N,p,si,taus,select.random,qcond1,
tau2par1,qcond2,tau2par2)
rimperfect.trivariateVineCopulaREMADA.beta(N,p,g,taus,select.random,qcond1,
tau2par1,qcond2,tau2par2)
}
\arguments{
\item{N}{sample size}
\item{p}{Vector \eqn{(\pi_{1},\pi_{2},\pi_{3},\pi_{4},\pi_{5})}, where \eqn{\pi_1}
 is the meta-analytic parameter for the prevalence, \eqn{\pi_2} and \eqn{\pi_3} are the meta-analytic parameters for the sensitivity of the index and the reference test, respectively, and \eqn{\pi_4} and \eqn{\pi_5} are the meta-analytic parameters for the specificity of the index and the reference test, respectively.}
\item{si}{Vector \eqn{(\sigma_{1},\sigma_{2},\sigma_{3})}, where  \eqn{\sigma_t,\,t=1,\ldots,3} denote the between-study heterogeneities (normal margins)}
\item{g}{Vector \eqn{(\gamma_{1},\gamma_{2},\gamma_{3})} where\eqn{\gamma_t,\,t=1,\ldots,3} denote the between-study heterogeneities (beta margins)}
\item{taus}{Kendall's tau values}
\item{select.random}{vector \eqn{(t_{1},t_{2},t_3)}, where \eqn{1\leq t_1<t_2<t_3\leq 5}}
\item{qcond1}{function for the inverse of conditional copula cdf for the \eqn{(t_{1},t_{2})} bivariate margin; choices are \code{qcondbvn}, \code{qcondfrk}, \code{qcondcln}, \code{qcondcln90}, \code{qcondcln180} and \code{qcondcln270}}
\item{tau2par1}{function for maping Kendall's tau to copula parameter for the \eqn{(t_{1},t_{2})} bivariate margin; choices are \code{tau2par.bvn}, \code{tau2par.frk}, \code{tau2par.cln}, \code{tau2par.cln90}, \code{tau2par.cln180} and \code{tau2par.cln270}}
\item{qcond2}{function for the inverse of conditional copula cdf for the \eqn{(t_{2},t_{3})} bivariate margin; choices are \code{qcondbvn}, \code{qcondfrk}, \code{qcondcln}, \code{qcondcln90}, \code{qcondcln180} and \code{qcondcln270}}
\item{tau2par2}{function for maping Kendall's tau to copula parameter for the \eqn{(t_{2},t_{3})} bivariate margin; choices are  \code{tau2par.bvn}, \code{tau2par.frk}, \code{tau2par.cln}, \code{tau2par.cln90}, \code{tau2par.cln180} and \code{tau2par.cln270}}
}

\value{Simulated data with 4 columns and \eqn{N} rows. \describe{
    \item{y11}{the number of the test results where the index test outcome is positive and the  reference test outcome is positive}
\item{y10}{the number of the test results where the index test outcome is positive and the  reference test outcome is negative}
\item{y01}{the number of the test results where the index test outcome is negative and the  reference test outcome is positive}
\item{y00}{the number of the test results where the index test outcome is negative and the  reference test outcome is negative}
}
}
\references{
Nikoloulopoulos, A.K. (2025) Vine copula mixed models for meta-analysis of diagnostic accuracy studies without a gold standard. \emph{Biometrics}, \bold{81}(2), ujaf037. \doi{10.1093/biomtc/ujaf037}.
}

\examples{
N=59
p=c(0.631,0.653,0.902,0.843,0.987)
si=c(1.513,1.341,1.341)
taus=c(0.3,-0.3)
select.random=c(1,2,4)

out=rimperfect.trivariateVineCopulaREMADA.norm(N,p,si,taus,select.random,
qcondcln180,tau2par.cln180,qcondcln270,tau2par.cln270)

g=c(0.290,0.244,0.190)
out=rimperfect.trivariateVineCopulaREMADA.beta(N,p,g,taus,select.random,
qcondcln180,tau2par.cln180,qcondcln270,tau2par.cln270)
}
\keyword{simulation}

