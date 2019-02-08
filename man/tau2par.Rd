\name{tau2par}
\Rdversion{1.1}
\alias{tau2par}
\alias{tau2par.bvn}
\alias{tau2par.frk}
\alias{tau2par.cln}
\alias{tau2par.cln90}
\alias{tau2par.cln180}
\alias{tau2par.cln270}

\title{
Mapping of Kendall's tau and copula parameter 
}
\description{
Bivariate copulas: mapping of Kendall's tau and copula parameter.
}
\usage{
tau2par.bvn(tau)
tau2par.frk(tau)
tau2par.cln(tau)
tau2par.cln90(tau)
tau2par.cln180(tau)
tau2par.cln270(tau)
}
\arguments{
\item{tau}{Kendall's tau for the copula }
}

\details{
For abbreviations of names of copula families (after the dot in function
names), see \code{\link{dcop}} help page.

}
\value{
copula parameter
}
\references{
Joe H (1997) \emph{Multivariate Models and Dependence Concepts}. Chapman & Hall

Joe H (2014) \emph{Dependence Modeling with Copulas}. Chapman & Hall/CRC.

Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. Software for book: \emph{Dependence Modeling with Copulas}, Chapman & Hall/CRC, 2014.
}
\seealso{
\code{\link{dcop}}
}

\keyword{copula}
\keyword{dependence}
