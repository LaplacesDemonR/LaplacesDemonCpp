#include "interval.h"

using namespace Rcpp;
using namespace arma;

SEXP interval(SEXP X, SEXP A, SEXP B, SEXP REFLECT) {
  Rcpp::NumericVector x(X);
  double a = as<double>(A);
  double b = as<double>(B);
  bool reflect = REFLECT;
  int N = x.size();
  int out = 1;
  if (a > b) a = b;
  if (reflect == false) {
    for (int i = 1; i < N; i++) {
      if (x[i] < a) x[i] = a;
      else if (x[i] > b) x[i] = b;
    }
  }
  else {
    for (int i = 1; i < N; i++) {
      out = 1;
      while (out == 1) {
        if (x[i] >= a && x[i] <= b) out = 0;
	else if (x[i] < a) x[i] = a + a - x[i];
        else x[i] = b + b - x[i];
      }
    }
  }
  return wrap(x);
}

SEXP intervala(SEXP X, SEXP DIM, SEXP A, SEXP B, SEXP REFLECT) {
  Rcpp::NumericVector x(X);
  Rcpp::IntegerVector dims(DIM);
  double a = as<double>(A);
  double b = as<double>(B);
  bool reflect = REFLECT;
  int N = x.size();
  int out = 1;
  if (a > b) a = b;
  if (reflect == false) {
    for (int i = 0; i < N; i++) {
      if (x[i] < a) x[i] = a;
      else if (x[i] > b) x[i] = b;
    }
  }
  else {
    for (int i = 0; i < N; i++) {
      out = 1;
      while (out == 1) {
        if (x[i] >= a && x[i] <= b) out = 0;
        else if (x[i] < a) x[i] = a + a - x[i];
        else x[i] = b + b - x[i];
      }
    }
  }
  x.attr("dim") = dims;
  return wrap(x);
}
