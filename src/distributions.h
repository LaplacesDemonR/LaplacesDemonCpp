#ifndef _LaplacesDemonCpp_dhalfcauchy_H
#define _LaplacesDemonCpp_dhalfcauchy_H

#include <RcppArmadillo.h>

RcppExport SEXP dhalfcauchy(SEXP X, SEXP SCALE, SEXP LOGD);

#endif

#ifndef _LaplacesDemonCpp_dinvwishart_H
#define _LaplacesDemonCpp_dinvwishart_H

#include <RcppArmadillo.h>

RcppExport SEXP dinvwishart(SEXP SIGMA, SEXP NU, SEXP s, SEXP LOGD);

#endif

#ifndef _LaplacesDemonCpp_dinvwishartc_H
#define _LaplacesDemonCpp_dinvwishartc_H

#include <RcppArmadillo.h>

RcppExport SEXP dinvwishartc(SEXP u, SEXP NU, SEXP s, SEXP LOGD);

#endif

#ifndef _LaplacesDemonCpp_dmvn_H
#define _LaplacesDemonCpp_dmvn_H

#include <RcppArmadillo.h>

RcppExport SEXP dmvn(SEXP X, SEXP MU, SEXP SIGMA, SEXP LOGD);

#endif

#ifndef _LaplacesDemonCpp_dmvnc_H
#define _LaplacesDemonCpp_dmvnc_H

#include <RcppArmadillo.h>

RcppExport SEXP dmvnc(SEXP X, SEXP MU, SEXP u, SEXP LOGD);

#endif

#ifndef _LaplacesDemonCpp_dmvnp_H
#define _LaplacesDemonCpp_dmvnp_H

#include <RcppArmadillo.h>

RcppExport SEXP dmvnp(SEXP X, SEXP MU, SEXP OMEGA, SEXP LOGD);

#endif

#ifndef _LaplacesDemonCpp_dmvnpc_H
#define _LaplacesDemonCpp_dmvnpc_H

#include <RcppArmadillo.h>

RcppExport SEXP dmvnpc(SEXP X, SEXP MU, SEXP u, SEXP LOGD);

#endif

#ifndef _LaplacesDemonCpp_dwishart_H
#define _LaplacesDemonCpp_dwishart_H

#include <RcppArmadillo.h>

RcppExport SEXP dwishart(SEXP OMEGA, SEXP NU, SEXP s, SEXP LOGD);

#endif

#ifndef _LaplacesDemonCpp_dwishartc_H
#define _LaplacesDemonCpp_dwishartc_H

#include <RcppArmadillo.h>

RcppExport SEXP dwishartc(SEXP u, SEXP NU, SEXP s, SEXP LOGD);

#endif

#ifndef _LaplacesDemonCpp_rhalfcauchy_H
#define _LaplacesDemonCpp_rhalfcauchy_H

#include <RcppArmadillo.h>

RcppExport SEXP rhalfcauchy(SEXP N, SEXP SCALE);

#endif

#ifndef _LaplacesDemonCpp_rmvn_H
#define _LaplacesDemonCpp_rmvn_H

#include <RcppArmadillo.h>

RcppExport SEXP rmvn(SEXP MU, SEXP SIGMA);

#endif

#ifndef _LaplacesDemonCpp_rmvnp_H
#define _LaplacesDemonCpp_rmvnp_H

#include <RcppArmadillo.h>

RcppExport SEXP rmvnp(SEXP MU, SEXP OMEGA);

#endif

#ifndef _LaplacesDemonCpp_rmvnc_H
#define _LaplacesDemonCpp_rmvnc_H

#include <RcppArmadillo.h>

RcppExport SEXP rmvnc(SEXP MU, SEXP u);

#endif

#ifndef _LaplacesDemonCpp_rmvnpc_H
#define _LaplacesDemonCpp_rmvnpc_H

#include <RcppArmadillo.h>

RcppExport SEXP rmvnpc(SEXP MU, SEXP u);

#endif

#ifndef _LaplacesDemonCpp_rwishart_H
#define _LaplacesDemonCpp_rwishart_H

#include <RcppArmadillo.h>

RcppExport SEXP rwishart(SEXP NU, SEXP s);

#endif

#ifndef _LaplacesDemonCpp_rwishartc_H
#define _LaplacesDemonCpp_rwishartc_H

#include <RcppArmadillo.h>

RcppExport SEXP rwishartc(SEXP NU, SEXP s);

#endif
