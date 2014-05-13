#ifndef _LaplacesDemonCpp_interval_H
#define _LaplacesDemonCpp_interval_H

#include <RcppArmadillo.h>

RcppExport SEXP interval(SEXP X, SEXP A, SEXP B, SEXP REFLECT);

#endif

#ifndef _LaplacesDemonCpp_intervala_H
#define _LaplacesDemonCpp_intervala_H

#include <RcppArmadillo.h>

RcppExport SEXP intervala(SEXP X, SEXP DIM, SEXP A, SEXP B, SEXP REFLECT);

#endif
