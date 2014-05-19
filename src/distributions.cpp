#include "distributions.h"

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

SEXP rhalfcauchy(SEXP N, SEXP SCALE) {
  int n = as<int>(N);
  double scale = as<double>(SCALE);
  RNGScope scope;
  Rcpp::NumericVector p(M_PI * runif(n) / 2.0), x(n);
  for (int i = 0; i < n; i++) {
    x[i] = scale * tan(p[i]);
  }
  return wrap(x);
}

/*------------------------------------------------------------------------/
/ Inverse Wishart Distribution                                            /
/------------------------------------------------------------------------*/

SEXP dinvwishart(SEXP SIGMA, SEXP NU, SEXP s, SEXP LOGD) {
  arma::mat Sigma = as<arma::mat>(SIGMA);
  double nu = as<double>(NU);
  arma::mat S = as<arma::mat>(s);
  bool logd = as<bool>(LOGD);
  int k = Sigma.n_rows;
  double dens = 0, gamsum = 0;
  for (int i = 0; i < k; i++) {
    gamsum += lgamma((nu + 1.0 - i - 1.0) / 2.0);
  }
  dens = -((nu * k) / 2.0) * log(2.0) - ((k * (k - 1.0)) / 4.0) * 
    log(M_PI) - gamsum + (nu / 2.0) * log(arma::det(S)) - 
    ((nu + k + 1.0) / 2.0) * log(arma::det(Sigma)) - 
    0.5 * arma::trace(S * arma::inv(Sigma));
  if (logd == false) dens = exp(dens);
  return wrap(dens);
}

/*------------------------------------------------------------------------/
/ Inverse Wishart Distribution (Cholesky Parameterization)                /
/------------------------------------------------------------------------*/

SEXP dinvwishartc(SEXP u, SEXP NU, SEXP s, SEXP LOGD) {
  arma::mat U = as<arma::mat>(u);
  double nu = as<double>(NU);
  arma::mat S = as<arma::mat>(s);
  bool logd = as<bool>(LOGD);
  int k = U.n_rows;
  double dens = 0, gamsum = 0;
  arma::mat Sigma = arma::trans(U) * U;
  for (int i = 0; i < k; i++) {
    gamsum += lgamma((nu + 1.0 - i - 1.0) / 2.0);
  }
  dens = -((nu * k) / 2.0) * log(2.0) - ((k * (k - 1.0)) / 4.0) * 
    log(M_PI) - gamsum + (nu / 2.0) * log(arma::det(S)) - 
    ((nu + k + 1.0) / 2.0) * log(arma::det(Sigma)) - 
    0.5 * arma::trace(S * arma::inv(Sigma));
  if (logd == false) dens = exp(dens);
  return wrap(dens);
}

/*------------------------------------------------------------------------/
/ Laplace Distribution                                                    /
/------------------------------------------------------------------------*/

SEXP dlaplace(SEXP X, SEXP LOCATION, SEXP SCALE, SEXP LOGD) {
  Rcpp::NumericVector x(X), location(LOCATION), scale(SCALE);
  int N = max(Rcpp::NumericVector::create(x.size(), location.size(), 
    scale.size()));
  Rcpp::NumericVector xn = rep_len(x, N), locationn = rep_len(location, N), 
    scalen = rep_len(scale, N), dens(N);
  bool logd = as<bool>(LOGD);
  Rcpp::NumericVector xloc = xn - locationn;
  xloc = abs(xloc);
  for (int i = 0; i < N; i++) {
    dens[i] = (-xloc[i] / scalen[i]) - log(2 * scalen[i]);
  }
  if (logd == false) dens = exp(dens);
  return wrap(dens);
}

SEXP rlaplace(SEXP N, SEXP LOCATION, SEXP SCALE) {
  int n = as<int>(N);
  double signr = 0;
  Rcpp::NumericVector location(LOCATION), scale(SCALE);
  Rcpp::NumericVector locationn = rep_len(location, n),
    scalen = rep_len(scale, n);
  Rcpp::NumericVector r = runif(n), r2(r), x(n);
  RNGScope scope;
  for (int i = 0; i < n; i++) {
    if (r[i] > 0.5) {
      r2[i] = 1 - r[i];
      signr = 1;
    } else signr = -1;
    x[i] = locationn[i] - signr * scalen[i] * log(2 * r2[i]);
  }
  return wrap(x);
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
  for (int i = 0; i < n; i++) {
    arma::vec z = rooti * arma::trans(x.row(i) - mu.row(i));
    dens[i] = constants + rootisum - 0.5 * arma::sum(z%z);
  }
  if (logd == false) dens = exp(dens);
  return(wrap(dens));
}

SEXP rmvn(SEXP MU, SEXP SIGMA) {
  arma::mat mu = as<arma::mat>(MU);
  arma::mat Sigma = as<arma::mat>(SIGMA);
  int n = mu.n_rows;
  RNGScope scope;
  arma::mat z = arma::randn(n, Sigma.n_cols);
  return wrap(mu + z * arma::chol(Sigma));
}

/*------------------------------------------------------------------------/
/ Multivariate Normal Distribution (Cholesky Parameterization)            /
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

SEXP rmvnc(SEXP MU, SEXP u) {
  arma::mat mu = as<arma::mat>(MU);
  arma::mat U = as<arma::mat>(u);
  int n = mu.n_rows;
  RNGScope scope;
  arma::mat z = arma::randn(n, U.n_cols);
  return wrap(mu + z * U);
}

/*------------------------------------------------------------------------/
/ Multivariate Normal Distribution (Precision Parameterization)           /
/------------------------------------------------------------------------*/

SEXP dmvnp(SEXP X, SEXP MU, SEXP OMEGA, SEXP LOGD) {
  arma::mat x = as<arma::mat>(X);
  arma::mat mu = as<arma::mat>(MU);
  arma::mat Omega = as<arma::mat>(OMEGA);
  bool logd = as<bool>(LOGD);
  int k = x.n_cols;
  int n = x.n_rows;
  double logdetOmega = log(arma::det(Omega));
  double log2pi = std::log(2.0 * M_PI);
  arma::mat ss = x - mu;
  arma::mat z = ss * Omega % ss;
  Rcpp::NumericVector dens(n);
  for (int i = 0; i < n; i++) {
    dens[i] = (-k / 2.0) * log2pi + 0.5 * logdetOmega - 
      0.5 * arma::sum(z.row(i));
  }
  if (logd == false) dens = exp(dens);
  return wrap(dens);
}

SEXP rmvnp(SEXP MU, SEXP OMEGA) {
  arma::mat mu = as<arma::mat>(MU);
  arma::mat Omega = as<arma::mat>(OMEGA);
  int n = mu.n_rows;
  RNGScope scope;
  arma::mat z = arma::randn(n, Omega.n_cols);
  return wrap(mu + z * arma::chol(arma::inv(Omega)));
}

/*------------------------------------------------------------------------/
/ Multivariate Normal Distribution (Precision-Cholesky Parameterization)  /
/------------------------------------------------------------------------*/

SEXP dmvnpc(SEXP X, SEXP MU, SEXP u, SEXP LOGD) {
  arma::mat x = as<arma::mat>(X);
  arma::mat mu = as<arma::mat>(MU);
  arma::mat U = as<arma::mat>(u);
  arma::mat Omega = arma::trans(U) * U;
  bool logd = as<bool>(LOGD);
  int k = x.n_cols;
  int n = x.n_rows;
  double logdetOmega = log(arma::det(Omega));
  double log2pi = std::log(2.0 * M_PI);
  arma::mat ss = x - mu;
  arma::mat z = ss * Omega % ss;
  Rcpp::NumericVector dens(n);
  for (int i = 0; i < n; i++) {
    dens[i] = (-k / 2.0) * log2pi + 0.5 * logdetOmega - 
      0.5 * arma::sum(z.row(i));
  }
  if (logd == false) dens = exp(dens);
  return wrap(dens);
}

SEXP rmvnpc(SEXP MU, SEXP u) {
  arma::mat mu = as<arma::mat>(MU);
  arma::mat U = as<arma::mat>(u);
  int n = mu.n_rows;
  RNGScope scope;
  arma::mat z = arma::randn(n, U.n_cols);
  return wrap(mu + z * arma::inv(arma::trans(U)));
}

/*------------------------------------------------------------------------/
/ Multivariate t Distribution                                             /
/------------------------------------------------------------------------*/

SEXP dmvt(SEXP X, SEXP MU, SEXP s, SEXP DF, SEXP LOGD) {
  arma::mat x = as<arma::mat>(X);
  arma::mat mu = as<arma::mat>(MU);
  arma::mat S = as<arma::mat>(s);
  double df = as<double>(DF);
  bool logd = as<bool>(LOGD);
  int k = x.n_cols;
  int n = x.n_rows;
  double logdetS = log(arma::det(S));
  arma::mat ss = x - mu;
  arma::mat Omega = arma::inv(S);
  arma::mat z = ss * Omega % ss;
  Rcpp::NumericVector dens(n);
  for (int i = 0; i < n; i++) {
    dens[i] = lgamma((df + k) / 2.0) - lgamma(df / 2.0) + (k / 2.0)*df +
      (k / 2.0) * log(M_PI) + 0.5 * logdetS + ((df + k) / 2.0) *
      log(1.0 + (1.0 / df) * arma::sum(z.row(i)));
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
    gamsum += lgamma((nu + 1.0 - i - 1.0) / 2.0);
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

/*------------------------------------------------------------------------/
/ Wishart Distribution (Cholesky Parameterization)                        /
/------------------------------------------------------------------------*/

SEXP dwishartc(SEXP u, SEXP NU, SEXP s, SEXP LOGD) {
  arma::mat U = as<arma::mat>(u);
  arma::mat Omega = arma::trans(U) * U;
  double nu = as<double>(NU);
  arma::mat S = as<arma::mat>(s);
  bool logd = as<bool>(LOGD);
  int k = Omega.n_rows;
  double dens = 0, gamsum = 0;
  for (int i = 0; i < k; i++) {
    gamsum += lgamma((nu + 1.0 - i - 1.0) / 2.0);
  }
  dens = -((nu * k) / 2.0) * log(2.0) - ((k * (k - 1.0)) / 4.0) * 
    log(M_PI) - gamsum - (nu / 2.0) * log(arma::det(S)) + 
    ((nu - k - 1.0) / 2.0) * log(arma::det(Omega)) - 
    (arma::trace(arma::inv(S) * Omega) / 2.0);
  if (logd == false) dens = exp(dens);
  return wrap(dens);
}

SEXP rwishartc(SEXP NU, SEXP s) {
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
