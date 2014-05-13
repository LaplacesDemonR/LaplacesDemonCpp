###########################################################################
# Matrices                                                                #
#                                                                         #
# These are utility functions for matrices.                               #
###########################################################################

is.positive.definite <- function(x) {
     return(.Call("is_positive_definite", x, PACKAGE="LaplacesDemonCpp"))}

is.symmetric.matrix <- function(x) {
     return(.Call("is_symmetric_matrix", x, PACKAGE="LaplacesDemonCpp"))}

tr <- function(x) {return(.Call("tr", x, PACKAGE="LaplacesDemonCpp"))}

#End
