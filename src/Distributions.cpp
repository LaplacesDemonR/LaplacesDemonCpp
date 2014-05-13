#include "Distributions.h"

using namespace Rcpp;
using namespace arma;

/*------------------------------------------------------------------------/
/ Half-Cauchy Distribution                                                /
/------------------------------------------------------------------------*/

SEXP dhalfcauchy(SEXP X, SEXP SCALE, SEXP LOGD) {
  Rcpp::NumericVector x(X), scale(SCALE);
  bool logd = LOGD;
  int N = x.size(), Ns = scale.size();
  double sigma = 0;
  Rcpp::NumericVector dens(N);
  for (int i = 0; i < N; i++) {
    if(Ns < N) sigma = scale[0];
    else sigma = scale[i];
    dens[i] = log(2 * sigma) - log(M_PI * (pow(x[i], 2) + pow(sigma, 2)));
  }
  if (logd == false) dens = exp(dens);
  return wrap(dens);
}

/*------------------------------------------------------------------------/
/ Wishart Distribution                                                    /
/------------------------------------------------------------------------*/

SEXP dwishart(SEXP OMEGA, SEXP NU, SEXP s, SEXP LOGD) {
  arma::mat Omega = as<arma::mat>(OMEGA);
  double nu = as<double>(NU);
  arma::mat S = as<arma::mat>(s);
  bool logd = LOGD;
  int k = Omega.n_rows;
  double dens = 0, gamsum = 0;
  for (int i = 0; i < k; i++) {
    gamsum += lgamma((nu + 1 - i) / 2.0);
  }
  dens = -((nu * k) / 2.0) * log(2.0) - ((k * (k - 1.0)) / 4.0) * 
    log(M_PI) - gamsum - (nu / 2.0) * log(arma::det(S)) + 
    ((nu - k - 1.0) / 2.0) * log(arma::det(Omega)) - 
    (arma::trace(arma::inv(S) * Omega) / 2.0);
  if (logd == false) dens = exp(dens);
  return wrap(dens);
}

SEXP rwishart(SEXP NU, SEXP s) {
  Rcpp::NumericVector nu(NU);
  arma::mat S = as<arma::mat>(s);
  arma::mat CC = arma::chol(S);
  int n = S.n_cols;
  arma::mat a = arma::randn(n, n);
  for (int i = 0; i < n; i++) {
    a.diag()[i] = sqrt(as<double>(rchisq(1, nu[i])));
  }
  a = arma::trimatu(a);
  a = a * CC;
  a = arma::trans(a) * a;
  return wrap(a);
}
