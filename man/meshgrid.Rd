\name{meshgrid}
\alias{meshgrid}
\title{
  Generate a Mesh Grid
}
\description{
  Generate two matrices for use in three-dimensional plots.
}
\usage{
meshgrid(x, y = x)
}
\arguments{
  \item{x}{numerical vector, represents points along the x-axis.}
  \item{y}{numerical vector, represents points along the y-axis.}
}
\details{
  The rows of the output array X are copies of the vector x;
  columns of the output array Y are copies of the vector y.
}
\value{
  Returns two matrices as a list with \code{X} and \code{Y} components.
}
\references{ Hans W. Borchers. \emph{pracma: Practical Numerical Math Functions} 2015. URL
http://CRAN.R-project.org/package=pracma, R package version 1.8.3.
}
\examples{
meshgrid(1:5)$X
meshgrid(c(1, 2, 3), c(11, 12))
}
\keyword{ array }
