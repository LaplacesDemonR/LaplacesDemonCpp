#include "Distributions.h"

using namespace Rcpp;
using namespace arma;

/*------------------------------------------------------------------------/
/ Half-Cauchy Distribution                                                /
/------------------------------------------------------------------------*/

SEXP dhalfcauchy(SEXP X, SEXP SCALE, SEXP LOGD) {
  Rcpp::NumericVector x(X), scale(SCALE);
  int N = max(Rcpp::NumericVector::create(x.size(), scale.size()));
  Rcpp::NumericVector xn = rep_len(x, N), scalen = rep_len(scale, N);
  bool logd = as<bool>(LOGD);
  bool scaletest = any(scale <= 0).is_true();
  if (scaletest == true) stop("The scale parameter must be positive");
  Rcpp::NumericVector dens(N);
  for (int i = 0; i < N; i++) {
    dens[i] = log(2 * scalen[i]) - log(M_PI * (pow(xn[i], 2) + 
      pow(scalen[i], 2)));
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
  bool logd = as<bool>(LOGD);
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
