#include "Math.h"

using namespace Rcpp;
using namespace arma;

SEXP partial(SEXP MODEL, SEXP PARM, SEXP DATA, SEXP INTERVAL) {
  Rcpp::Function Model = MODEL;
  Rcpp::NumericVector parm(PARM);
  Rcpp::List Data = DATA;
  double Interval = as<double>(INTERVAL);
  int N = parm.size(), fin = 0;
  double partial = 0;
  Rcpp::NumericVector df(N);
  Rcpp::List Mold = Model(parm, Data);
  Rcpp::List Mnew = clone(Mold);
  for (int i = 0; i < N; i++) {
    Rcpp::List Mold2 = clone(Mold);
    Rcpp::NumericVector dx = Mold2["parm"];
    dx[i] += Interval;
    Mnew = Model(dx, Data);
    partial = (as<double>(Mnew["LP"]) - as<double>(Mold["LP"])) / Interval;
    fin = ::R_finite(partial);
    if (fin == 1) df[i] = partial;
    else df[i] = 0.0;
  }
  return wrap(df);
}
