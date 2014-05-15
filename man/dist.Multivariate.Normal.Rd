\name{dist.Multivariate.Normal}
\alias{dmvn}
\alias{rmvn}
\title{Multivariate Normal Distribution}
\description{
  These are C++ versions of the \code{dmvn} and \code{rmvn} functions in
  the \pkg{LaplacesDemon} package.
}
\usage{
dmvn(x, mu, Sigma, log=FALSE)
rmvn(n=1, mu=rep(0,ncol(Sigma)), Sigma)
}
\arguments{
  \item{x}{See \code{dmvn}.}
  \item{n}{See \code{rmvn}.}
  \item{mu}{See \code{dmvn}.}
  \item{Sigma}{See \code{dmvn}.}
  \item{log}{See \code{dmvn}.}
}
\details{
  See \code{dmvn}.
}
\value{
  See \code{dmvn}.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \url{www.bayesian-inference.com/softwaredoc/dist.Multivariate.Normal}
}
\keyword{Distribution}