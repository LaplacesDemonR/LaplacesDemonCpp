###########################################################################
# interval                                                                #
#                                                                         #
# The purpose of the interval function is to constrain the element(s) of  #
# a scalar, vector, matrix, or array to the interval [a,b].               #
###########################################################################

interval <- function(x, a=-Inf, b=Inf, reflect=TRUE)
     {
     if(missing(x)) stop("The x argument is required.")
     if(a > b) stop("a > b.")
     if(is.vector(x) == TRUE)
          return(.Call("interval", x, a, b, reflect,
               PACKAGE="LaplacesDemonCpp"))
     else if(is.array(x) == TRUE)
          return(.Call("intervala", x, dim(x), a, b, reflect,
               PACKAGE="LaplacesDemonCpp"))
     } 

#End
