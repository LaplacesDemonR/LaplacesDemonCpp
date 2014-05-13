#ifndef _LaplacesDemonCpp_dhalfcauchy_H
#define _LaplacesDemonCpp_dhalfcauchy_H

#include <RcppArmadillo.h>

RcppExport SEXP dhalfcauchy(SEXP X, SEXP SCALE, SEXP LOGD);

#endif

#ifndef _LaplacesDemonCpp_dwishart_H
#define _LaplacesDemonCpp_dwishart_H

#include <RcppArmadillo.h>

RcppExport SEXP dwishart(SEXP OMEGA, SEXP NU, SEXP s, SEXP LOGD);

#endif


#ifndef _LaplacesDemonCpp_rwishart_H
#define _LaplacesDemonCpp_rwishart_H

#include <RcppArmadillo.h>

RcppExport SEXP rwishart(SEXP NU, SEXP s);

#endif
