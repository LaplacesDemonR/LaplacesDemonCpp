###########################################################################
# Half-Cauchy Distribution                                                #
###########################################################################

dhalfcauchy <- function(x, scale=25, log=FALSE)
     {
     return(.Call("dhalfcauchy", x, scale, log,
          PACKAGE="LaplacesDemonCpp"))
     }

###########################################################################
# Wishart Distribution                                                    #
###########################################################################

dwishart <- function(Omega, nu, S, log=FALSE)
     {
     return(.Call("dwishart", Omega, nu, S, log,
          PACKAGE="LaplacesDemonCpp"))
     }

rwishart <- function(nu, S)
     {
     return(.Call("rwishart", nu, S, PACKAGE="LaplacesDemonCpp"))
     }
