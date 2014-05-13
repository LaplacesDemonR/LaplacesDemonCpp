###########################################################################
# Matrices                                                                #
#                                                                         #
# These are utility functions for matrices.                               #
###########################################################################

is.symmetric.matrix <- function(x) {
     return(.Call("is_symmetric_matrix", x, PACKAGE="LaplacesDemonCpp"))}

tr <- function(x) {return(.Call("tr", x, PACKAGE="LaplacesDemonCpp"))}

#End