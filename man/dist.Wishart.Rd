\name{dist.Wishart}
\alias{dwishart}
\alias{rwishart}
\title{Wishart Distribution}
\description{
  This is a C++ version of the \code{dwishart} and \code{rwishart}
  functions in the \pkg{LaplacesDemon} package.
}
\usage{
   dwishart(Omega, nu, S, log=FALSE)
   rwishart(nu, S)
}
\arguments{
  \item{Omega}{See \code{dwishart}.}
  \item{nu}{See \code{dwishart}.}
  \item{S}{See \code{dwishart}.}
  \item{log}{See \code{dwishart}.}
}
\details{
  See \code{dwishart}.

  Results do not exactly match \pkg{LaplacesDemon}. Bug? Please send an
  email if you solve it.
}
\value{
  See \code{dwishart}.
}
\seealso{
  \url{www.bayesian-inference.com/softwaredoc/dist.Wishart}
}
\keyword{Distribution}
