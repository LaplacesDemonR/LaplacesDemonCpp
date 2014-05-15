###########################################################################
# Math                                                                    #
#                                                                         #
# This is a collection of functions to facilitate math.                   #
###########################################################################

partial <- function(Model, parm, Data, Interval=1e-6, Method="simple")
     {
     f <- Model(parm, Data)[["LP"]]
     n <- length(parm)
     if(Method == "simple") {
          return(.Call("partial", Model, parm, Data, Interval=1e-6,
               PACKAGE="LaplacesDemonCpp"))
          }
     else if(Method == "Richardson") {
          zero.tol <- sqrt(.Machine$double.eps / 7e-7)
          d <- 0.0001
          r <- 4
          v <- 2
          a <- matrix(NA, r, n)
          h <- abs(d*parm) + Interval*{abs(parm) < zero.tol}
          for (k in 1:r) {
               if(n == 1)
                    a[k,] <- {Model(parm + h, Data)[["LP"]] -
                         Model(parm - h, Data)[["LP"]]} / (2*h)
               else for (i in 1:n) {
                    if((k != 1) && {abs(a[(k-1),i]) < 1e-20}) a[k,i] <- 0
	            else a[k,i] <- (Model(parm + h*(i == seq(n)), Data)[["LP"]] - 
	                 Model(parm - h*(i == seq(n)), Data)[["LP"]]) / (2*h[i])}
               a[k,which(!is.finite(a[k,]))] <- 0
               h <- h / v}
          for (m in 1:(r - 1))
               a <- (a[2:(r+1-m),,drop=FALSE]*(4^m)-a[1:(r-m),,drop=FALSE])/(4^m-1)
          return(c(a))
          }
     else stop("The", Method, "method is unknown.")
     }

#End
