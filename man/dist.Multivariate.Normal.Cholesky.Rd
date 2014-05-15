\name{dist.Multivariate.Normal.Cholesky}
\alias{dmvnc}
\alias{rmvnc}
\title{Multivariate Normal Distribution: Cholesky Parameterization}
\description{
  These are C++ versions of the \code{dmvnc} and \code{rmvnc} functions
  in the \pkg{LaplacesDemon} package.
}
\usage{
dmvnc(x, mu, U, log=FALSE)
rmvnc(n=1, mu=rep(0,ncol(U)), U)
}
\arguments{
  \item{x}{See \code{dmvnc}.}
  \item{n}{See \code{rmvnc}.}
  \item{mu}{See \code{dmvnc}.}
  \item{U}{See \code{dmvnc}.}
  \item{log}{See \code{dmvnc}.}
}
\details{
  See \code{dmvnc}.
}
\value{
  See \code{dmvnc}.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \url{www.bayesian-inference.com/softwaredoc/dist.Multivariate.Normal.Cholesky}
}
\keyword{Distribution}