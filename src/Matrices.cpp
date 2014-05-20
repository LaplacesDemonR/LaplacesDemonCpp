#include "Matrices.h"

using namespace Rcpp;
using namespace arma;

SEXP covar(SEXP x) {
  arma:mat X = as<arma::mat>(x);
  int n =X.n_rows;
  arma::vec ones(n);
  arma::mat x = X - ones * arma::trans(ones) * (X/n);
  arma::mat covar = arma::trans(x) * x / (n-1);
  return wrap(covar);
}

SEXP is_positive_definite(SEXP x) {
  arma::mat X = as<arma::mat>(x);
  arma::vec eigs = arma::eig_sym(X);
  bool pd = false;
  if (all(eigs > 0)) pd = true;
  return wrap(pd);
  // Thankfully, this does return false with a complex number.
}

SEXP is_positive_semidefinite(SEXP x) {
  arma::mat X = as<arma::mat>(x);
  arma::vec eigs = arma::eig_sym(X);
  bool pd = false;
  if (all(eigs >= 0)) pd = true;
  return wrap(pd);
  // Thankfully, this does return false with a complex number.
}

SEXP is_symmetric_matrix(SEXP x) {
  arma::mat X = as<arma::mat>(x);
  arma::mat Xt = arma::trans(X);
  int same = 0, n = X.n_rows, n2 = pow(n, 2);
  bool symmetric = false;
  for (int i = 0; i < n; i++) {
    same += arma::sum(X.row(i) == Xt.row(i));
  }
  if (same == n2) symmetric = true;
  return wrap(symmetric);
}

SEXP tr(SEXP x) {
  arma::mat X = as<arma::mat>(x);
  double tra = arma::trace(X);
  return wrap(tra);
}
