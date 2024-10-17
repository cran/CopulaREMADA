\name{mutinom6dVineCopulaREMADA}
\Rdversion{1.1}
\alias{multinom6dVineCopulaREMADA}
\alias{multinom6dVineCopulaREMADA.norm}
\alias{multinom6dVineCopulaREMADA.beta}


\title{
Maximum likelhood estimation  for multinomial six-variate 1-truncated D-vine copula mixed models for meta-analysis of two diagnostic tests accounting for within and between studies dependence
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
multinom6dVineCopulaREMADA.norm(y001,y011,y101,y111,y000,y010,y100,y110,
          gl,mgrid,qcond1,qcond2,qcond3,qcond4,qcond5,
          tau2par1,tau2par2,tau2par3,tau2par4,tau2par5,
          sel1,sel2,sel3)
multinom6dVineCopulaREMADA.beta(y001,y011,y101,y111,y000,y010,y100,y110,
          gl,mgrid,qcond1,qcond2,qcond3,qcond4,qcond5,
          tau2par1,tau2par2,tau2par3,tau2par4,tau2par5,
          sel1,sel2,sel3)
}
\arguments{
\item{y001}{the number of the test results in the diseased where the test 1 outcome is negative and  the test 2 outcome is negative}
\item{y011}{the number of the test results in the diseased where the test 1 outcome is negative and  the test 2 outcome is positive}
\item{y101}{the number of the test results in the diseased where the test 1 outcome is positive and  the test 2 outcome is negative}
\item{y111}{the number of the test results in the diseased where the test 1 outcome is positive and  the test 2 outcome is positive}
\item{y000}{the number of the test results in the non-diseased where the test 1 outcome is negative and  the test 2 outcome is negative}
\item{y010}{the number of the test results in the non-diseased where the test 1 outcome is negative and  the test 2 outcome is positive}
\item{y100}{the number of the test results in the non-diseased where the test 1 outcome is positive and  the test 2 outcome is negative}
\item{y110}{the number of the test results in the non-diseased where the test 1 outcome is positive and  the test 2 outcome is positive}
\item{gl}{a list containing the components of Gauss-Legendre nodes \code{gl$nodes} and weights \code{gl$weights}}
\item{mgrid}{
a list containing six-dimensional arrays.  Replicates of the quadrature points  that produce a 6-dimensional full grid}
\item{qcond1}{function for the inverse conditional copula cdf at the (1,2) bivariate margin}
\item{qcond2}{function for the inverse conditional copula cdf at the (2,3) bivariate margin}
\item{qcond3}{function for the inverse conditional copula cdf at the (3,4) bivariate margin}
\item{qcond4}{function for the inverse conditional copula cdf at the (4,5) bivariate margin}
\item{qcond5}{function for the inverse conditional copula cdf at the (5,6) bivariate margin}
\item{tau2par1}{function for maping Kendall's tau at the (1,2) bivariate margin to copula parameter}
\item{tau2par2}{function for maping Kendall's tau at the (2,3) bivariate margin to copula parameter}
\item{tau2par3}{function for maping Kendall's tau at the (3,4) bivariate margin to copula parameter}
\item{tau2par4}{function for maping Kendall's tau at the (4,5) bivariate margin to copula parameter}
\item{tau2par5}{function for maping Kendall's tau at the (5,6) bivariate margin to copula parameter}
\item{sel1}{Indicates the position of bivariate copulas with positive dependence only such as the Clayton  and the Clayton rotated by 180 degrees}
\item{sel2}{Indicates the position of bivariate copulas with negative dependence only such as the Clayton rotated by 90 degrees and the Clayton rotated by 270 degrees}
\item{sel3}{Indicates the position of bivariate copulas with comprehensive dependence such as the BVN and Frank copulas}
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
Nikoloulopoulos, A.K. (2024) Joint meta-analysis  of two diagnostic tests accounting for within and between studies dependence. Statistical Methods in Medical Research. \doi{10.1177/09622802241269645}

}
\seealso{
\code{\link{rmultinom6dVineCopulaREMADA}}
}
\examples{
\donttest{
data(Down)
attach(Down)
y111=down_n_11 
y110=nodown_n_11
y101=down_n_10
y100=nodown_n_10 
y001=down_n_00 
y000=nodown_n_00 
y010=nodown_n_01
y011=down_n_01

nq=15
gl=gauss.quad.prob(nq,"uniform")
data(mgrid6d)

tau2par1=tau2par.cln90
qcond1=qcondcln90
tau2par3=tau2par4=tau2par5=tau2par.cln
qcond3=qcond4=qcond5=qcondcln
tau2par2=tau2par.bvn
qcond2=qcondbvn

sel1=3:5; sel2=1; sel3=2

est=multinom6dVineCopulaREMADA.norm(y001,y011,y101,y111,
y000,y010,y100,y110,gl,mgrid,qcond1,qcond2,qcond3,qcond4,qcond5,
tau2par1,tau2par2,tau2par3,tau2par4,tau2par5,sel1,sel2,sel3)

detach(Down)
}}
\keyword{maximum likelihood}
\keyword{copula}
\keyword{distribution}

