#ifndef _LaplacesDemonCpp_is_matrices_H
#define _LaplacesDemonCpp_is_matrices_H

#include <RcppArmadillo.h>

RcppExport SEXP is_positive_definite(SEXP x);

RcppExport SEXP is_positive_semidefinite(SEXP x);

RcppExport SEXP is_symmetric_matrix(SEXP x);

RcppExport SEXP tr(SEXP x);

#endif
