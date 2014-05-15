###########################################################################
# Half-Cauchy Distribution                                                #
###########################################################################

dhalfcauchy <- function(x, scale=25, log=FALSE)
     {
     if(any(scale <= 0)) stop("The scale parameter must be positive.")
     return(.Call("dhalfcauchy", as.vector(x), as.vector(scale), log,
          PACKAGE="LaplacesDemonCpp"))
     }

###########################################################################
# Multivariate Normal Distribution                                        #
###########################################################################

dmvn <- function(x, mu, Sigma, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rbind(mu)
     nmax <- max(nrow(x), nrow(mu))
     x <- x[rep(1:nrow(x), len=nmax),]
     mu <- mu[rep(1:nrow(mu), len=nmax),]
     if(missing(Sigma)) Sigma <- diag(ncol(x))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     Sigma <- as.symmetric.matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     return(.Call("dmvn", x, mu, Sigma, log, PACKAGE="LaplacesDemonCpp"))
     }

rmvn <- function(n=1, mu=rep(0,ncol(Sigma)), Sigma)
     {
     mu <- rbind(mu)
     if(n > nrow(mu)) mu <- mu[rep(1:nrow(mu), len=n),]
     if(missing(Sigma)) Sigma <- diag(ncol(mu))
     if(!is.matrix(Sigma)) Sigma <- matrix(Sigma)
     if(!is.positive.definite(Sigma))
          stop("Matrix Sigma is not positive-definite.")
     return(.Call("rmvn", mu, Sigma, PACKAGE="LaplacesDemonCpp"))
     }

###########################################################################
# Multivariate Normal Distribution (Cholesky Parameterization)            #
###########################################################################

dmvnc <- function(x, mu, U, log=FALSE)
     {
     if(!is.matrix(x)) x <- rbind(x)
     if(!is.matrix(mu)) mu <- rbind(mu)
     nmax <- max(nrow(x), nrow(mu))
     x <- x[rep(1:nrow(x), len=nmax),]
     mu <- mu[rep(1:nrow(mu), len=nmax),]
     if(missing(U)) stop("Upper triangular U is required.")
     return(.Call("dmvnc", x, mu, U, log, PACKAGE="LaplacesDemonCpp"))
     }

rmvnc <- function(n=1, mu=rep(0,ncol(U)), U)
     {
     mu <- rbind(mu)
     if(n > nrow(mu)) mu <- mu[rep(1:nrow(mu), len=n),]
     if(missing(U)) U <- diag(ncol(mu))
     if(!is.matrix(U)) U <- matrix(U)
     return(.Call("rmvnc", mu, U, PACKAGE="LaplacesDemonCpp"))
     }

###########################################################################
# Multivariate Normal Distribution (Precision Parameterization)           #
###########################################################################

rmvnp <- function(n=1, mu=rep(0,ncol(Omega)), Omega)
     {
     mu <- rbind(mu)
     if(n > nrow(mu)) mu <- mu[rep(1:nrow(mu), len=n),]
     if(missing(Omega)) Omega <- diag(ncol(mu))
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega))
          stop("Matrix Omega is not positive-definite.")
     return(.Call("rmvnp", mu, Omega, PACKAGE="LaplacesDemonCpp"))
     }

###########################################################################
# Multivariate Normal Distribution (Precision-Cholesky Parameterization)  #
###########################################################################

rmvnpc <- function(n=1, mu=rep(0,ncol(U)), U)
     {
     mu <- rbind(mu)
     if(n > nrow(mu)) mu <- mu[rep(1:nrow(mu), len=n),]
     if(missing(U)) U <- diag(ncol(mu))
     if(!is.matrix(U)) U <- matrix(U)
     return(.Call("rmvnpc", mu, U, PACKAGE="LaplacesDemonCpp"))
     }

###########################################################################
# Wishart Distribution                                                    #
###########################################################################

dwishart <- function(Omega, nu, S, log=FALSE)
     {
     if(!is.matrix(Omega)) Omega <- matrix(Omega)
     if(!is.positive.definite(Omega)) 
          stop("Matrix Omega is not positive-definite.")
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S))
          stop("Matrix S is not positive-semidefinite.")
     if(!identical(dim(Omega), dim(S))) 
          stop("The dimensions of Omega and S differ.")
     if(nu < nrow(S)) 
          stop("The nu parameter is less than the dimension of S.")
     return(.Call("dwishart", Omega, nu, S, log,
          PACKAGE="LaplacesDemonCpp"))
     }

rwishart <- function(nu, S)
     {
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S)) 
          stop("Matrix S is not positive-semidefinite.")
     if(nu < nrow(S))
          stop("The nu parameter is less than the dimension of S.")
     return(.Call("rwishart", nu, S, PACKAGE="LaplacesDemonCpp"))
     }

###########################################################################
# Wishart Distribution (Cholesky Parameterization)                        #
###########################################################################

dwishartc <- function(U, nu, S, log=FALSE)
     {
     if(!is.matrix(U)) U <- matrix(U)
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S))
          stop("Matrix S is not positive-semidefinite.")
     if(!identical(dim(U), dim(S))) 
          stop("The dimensions of U and S differ.")
     if(nu < nrow(S)) 
          stop("The nu parameter is less than the dimension of S.")
     return(.Call("dwishartc", U, nu, S, log, PACKAGE="LaplacesDemonCpp"))
     }

rwishartc <- function(nu, S)
     {
     if(!is.matrix(S)) S <- matrix(S)
     if(!is.positive.semidefinite(S)) 
          stop("Matrix S is not positive-semidefinite.")
     if(nu < nrow(S))
          stop("The nu parameter is less than the dimension of S.")
     return(.Call("rwishartc", nu, S, PACKAGE="LaplacesDemonCpp"))
     }

#End
