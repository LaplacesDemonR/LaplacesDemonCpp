#ifndef _LaplacesDemonCpp_distributions_H
#define _LaplacesDemonCpp_distributions_H

#include <RcppArmadillo.h>

/*------------------------------------------------------------------------/
/ Categorical Distribution                                                /
/------------------------------------------------------------------------*/

RcppExport SEXP rcat(SEXP P);

/*------------------------------------------------------------------------/
/ Half-Cauchy Distribution                                                /
/------------------------------------------------------------------------*/

RcppExport SEXP dhalfcauchy(SEXP X, SEXP SCALE, SEXP LOGD);

RcppExport SEXP rhalfcauchy(SEXP N, SEXP SCALE);

/*------------------------------------------------------------------------/
/ Inverse Wishart Distribution                                            /
/------------------------------------------------------------------------*/

RcppExport SEXP dinvwishart(SEXP SIGMA, SEXP NU, SEXP s, SEXP LOGD);

/*------------------------------------------------------------------------/
/ Inverse Wishart Distribution (Cholesky Parameterization)                /
/------------------------------------------------------------------------*/

RcppExport SEXP dinvwishartc(SEXP u, SEXP NU, SEXP s, SEXP LOGD);

/*------------------------------------------------------------------------/
/ Laplace Distribution                                                    /
/------------------------------------------------------------------------*/

RcppExport SEXP dlaplace(SEXP X, SEXP LOCATION, SEXP SCALE, SEXP LOGD);

RcppExport SEXP rlaplace(SEXP X, SEXP LOCATION, SEXP SCALE);

/*------------------------------------------------------------------------/
/ Multivariate Normal Distribution                                        /
/------------------------------------------------------------------------*/

RcppExport SEXP dmvn(SEXP X, SEXP MU, SEXP SIGMA, SEXP LOGD);

RcppExport SEXP rmvn(SEXP MU, SEXP SIGMA);

/*------------------------------------------------------------------------/
/ Multivariate Normal Distribution (Cholesky Parameterization)            /
/------------------------------------------------------------------------*/

RcppExport SEXP dmvnc(SEXP X, SEXP MU, SEXP u, SEXP LOGD);

RcppExport SEXP rmvnc(SEXP MU, SEXP u);

/*------------------------------------------------------------------------/
/ Multivariate Normal Distribution (Precision Parameterization)           /
/------------------------------------------------------------------------*/

RcppExport SEXP dmvnp(SEXP X, SEXP MU, SEXP OMEGA, SEXP LOGD);

RcppExport SEXP rmvnp(SEXP MU, SEXP OMEGA);

/*------------------------------------------------------------------------/
/ Multivariate Normal Distribution (Precision-Cholesky Parameterization)  /
/------------------------------------------------------------------------*/

RcppExport SEXP dmvnpc(SEXP X, SEXP MU, SEXP u, SEXP LOGD);

RcppExport SEXP rmvnpc(SEXP MU, SEXP u);

/*------------------------------------------------------------------------/
/ Multivariate t Distribution                                             /
/------------------------------------------------------------------------*/

RcppExport SEXP dmvt(SEXP X, SEXP MU, SEXP s, SEXP DF, SEXP LOGD);

/*------------------------------------------------------------------------/
/ Wishart Distribution                                                    /
/------------------------------------------------------------------------*/

RcppExport SEXP dwishart(SEXP OMEGA, SEXP NU, SEXP s, SEXP LOGD);

RcppExport SEXP rwishart(SEXP NU, SEXP s);

/*------------------------------------------------------------------------/
/ Wishart Distribution (Cholesky Parameterization)                        /
/------------------------------------------------------------------------*/

RcppExport SEXP dwishartc(SEXP u, SEXP NU, SEXP s, SEXP LOGD);

RcppExport SEXP rwishartc(SEXP NU, SEXP s);

#endif
