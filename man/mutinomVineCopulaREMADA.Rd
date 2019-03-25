\name{mutinomVineCopulaREMADA}
\Rdversion{1.1}
\alias{multinomVineCopulaREMADA}
\alias{multinomVineCopulaREMADA.norm}
\alias{multinomVineCopulaREMADA.beta}
\alias{tmultinomVineCopulaREMADA}
\alias{tmultinomVineCopulaREMADA.norm}
\alias{tmultinomVineCopulaREMADA.beta}

\title{
Maximum likelhood estimation  for multinomial quadrivariate (truncated) D-vine copula mixed models for diagnostic test accurracy studies accounting for non-evaluable outcomes
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
multinomVineCopulaREMADA.norm(TP,FN,FP,TN,NEP,NEN,
                                 gl,mgrid,qcond1,pcond1,tau2par1,
                                 qcond2,pcond2,tau2par2)
multinomVineCopulaREMADA.beta(TP,FN,FP,TN,NEP,NEN,
                                 gl,mgrid,qcond1,pcond1,tau2par1,
                                 qcond2,pcond2,tau2par2)
tmultinomVineCopulaREMADA.norm(TP,FN,FP,TN,NEP,NEN,
                                 gl,mgrid,
                                 qcond1,pcond1,tau2par1,
                                 qcond2,pcond2,tau2par2)
tmultinomVineCopulaREMADA.beta(TP,FN,FP,TN,NEP,NEN,
                                 gl,mgrid,
                                 qcond1,pcond1,tau2par1,
                                 qcond2,pcond2,tau2par2)
}
\arguments{
\item{TP}{the number of true positives}
\item{FN}{the number of false negatives}
\item{FP}{the number of false positives}
\item{TN}{the number of true negatives}
\item{NEP}{the number of non-evaluable positives in the presence of non-evaluable subjects}
\item{NEN}{the number of non-evaluable negatives in the presence of non-evaluable subjects}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{
a list containing 4-dimensional arrays.}
\item{qcond1}{function for the inverse conditional copula cdf at the (1,2) and (3,4) bivariate margin}
\item{pcond1}{function for the conditional copula cdf at the (1,2) and (3,4) bivariate margin}
\item{tau2par1}{function for maping Kendall's tau at the (1,2) and (3,4) bivariate margin to copula parameter}
\item{qcond2}{function for the inverse conditional copula cdf at the (2,3) bivariate margin}
\item{pcond2}{function for the conditional copula cdf at the (2,3) bivariate margin}
\item{tau2par2}{function for maping Kendall's tau at the (2,3) bivariate margin to copula parameter}
}
\value{A list containing the following components:

\item{minimum}{the value of the estimated minimum of the negative log-likelihood}
\item{estimate}{the MLE}
\item{gradient}{the gradient at the estimated minimum of of the negative log-likelihood}
\item{hessian}{the hessian at the estimated minimum of the negative log-likelihood}
\item{code}{an integer indicating why the optimization process terminated}
\item{iterations}{the number of iterations performed}

For more details see  \code{\link{nlm}}
}
\references{
Nikoloulopoulos, A.K. (2018) A multinomial quadrivariate D-vine copula mixed model for diagnostic studies meta-analysis accounting for non-evaluable subjects.  \emph{ArXiv e-prints}, arXiv:1812.05915. \url{https://arxiv.org/abs/1812.05915}.

}
\seealso{
\code{\link{rmultinomVineCopulaREMADA}}
}
\examples{
nq=15
gl=gauss.quad.prob(nq,"uniform")
data(mgrid15)

data(MK2016)
attach(MK2016)

#out=tmultinomVineCopulaREMADA.beta(TP,FN,FP,TN,NEP,NEN,
#gl,mgrid15,qcondcln180,pcondcln180,tau2par.cln180,
#qcondcln90,pcondcln90,tau2par.cln90)

detach(MK2016)
}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

