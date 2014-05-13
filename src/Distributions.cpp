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
  Rcpp::NumericVector dens(N);
  for (int i = 0; i < N; i++) {
    dens[i] = log(2 * scalen[i]) - log(M_PI * (pow(xn[i], 2) + 
      pow(scalen[i], 2)));
  }
  if (logd == false) dens = exp(dens);
  return wrap(dens);
}

/*------------------------------------------------------------------------/
/ Multivariate Normal Distribution                                        /
/------------------------------------------------------------------------*/

SEXP dmvn(SEXP X, SEXP MU, SEXP SIGMA, SEXP LOGD) {
  arma::mat x = as<arma::mat>(X);
  arma::mat mu = as<arma::mat>(MU);
  arma::mat Sigma = as<arma::mat>(SIGMA);
  bool logd = as<bool>(LOGD);
  int n = x.n_rows;
  int k = x.n_cols;
  Rcpp::NumericVector dens(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double log2pi = std::log(2.0 * M_PI);
  double constants = -(static_cast<double>(k)/2.0) * log2pi;
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans(x.row(i) - mu.row(i));
    dens[i] = constants + rootisum - 0.5 * arma::sum(z%z);
  }
  if (logd == false) dens = exp(dens);
  return(wrap(dens));
}

/*------------------------------------------------------------------------/
/ Multivariate Normal Distribution: Cholesky Parameterization             /
/------------------------------------------------------------------------*/

SEXP dmvnc(SEXP X, SEXP MU, SEXP u, SEXP LOGD) {
  arma::mat x = as<arma::mat>(X);
  arma::mat mu = as<arma::mat>(MU);
  arma::mat U = as<arma::mat>(u);
  bool logd = as<bool>(LOGD);
  int n = x.n_rows;
  int k = x.n_cols;
  Rcpp::NumericVector dens(n);
  arma::mat rooti = arma::trans(arma::inv(U));
  double rootisum = arma::sum(log(rooti.diag()));
  double log2pi = std::log(2.0 * M_PI);
  double constants = -(static_cast<double>(k)/2.0) * log2pi;
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans(x.row(i) - mu.row(i));
    dens[i] = constants + rootisum - 0.5 * arma::sum(z%z);
  }
  if (logd == false) dens = exp(dens);
  return(wrap(dens));
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
  arma::mat x = arma::randn(n, n);
  for (int i = 0; i < n; i++) {
    x.diag()[i] = sqrt(as<double>(rchisq(1, nu[i])));
  }
  x = arma::trimatu(x);
  x = x * CC;
  x = arma::trans(x) * x;
  return wrap(x);
}
