#include "LaplacesDemon.h"
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

SEXP mcmcamwg(SEXP MODEL, SEXP DATA, SEXP ITERATIONS, SEXP STATUS, 
  SEXP THINNING, SEXP SPECS, SEXP ACCEPTANCE, SEXP DEV, SEXP DIAGCOVAR, 
  SEXP liv, SEXP MON, SEXP MO0, SEXP SCALEF, SEXP THINNED, SEXP TUNING) {
  Function Model(MODEL);
  Rcpp::List Data(DATA), Specs(SPECS), Mo0(MO0);
  int Iterations = as<int>(ITERATIONS), Status = as<int>(STATUS),
    Thinning = as<int>(THINNING), LIV = as<int>(liv), 
    ScaleF = as<int>(SCALEF), Periodicity = as<int>(Specs["Periodicity"]),
    ACCEPT = as<int>(ACCEPTANCE);
  Rcpp::IntegerVector Acceptance(LIV, ACCEPT);
  Rcpp::NumericMatrix Dev(DEV), Mon(MON), thinned(THINNED);
  Rcpp::NumericVector tuning(TUNING);
  int a_iter = 0, t_iter = 0, fins = 0, mcols = Mon.ncol();
  double log_alpha = 0;
  Rcpp::NumericMatrix DiagCovar(floor(Iterations / Periodicity), LIV);
  Rcpp::List Mo1 = clone(Mo0);
  RNGScope scope;
  for (int iter = 0; iter < Iterations; iter++) {
    // Print Status
    if ((iter + 1) % Status == 0) {
      Rcpp::Rcout << "Iteration: " << iter + 1 << 
           ",   Proposal: Componentwise,   LP: " << 
           floor(as<double>(Mo0["LP"]) * 100) / 100 << std::endl;
    }
    // Save Thinned Samples
    if ((iter + 1) % Thinning == 0) {
      t_iter = floor((iter) / Thinning) + 1;
      thinned(t_iter,_) = as<Rcpp::NumericVector>(Mo0["parm"]);
      Dev(t_iter,_) = as<Rcpp::NumericVector>(Mo0["Dev"]);
      Mon(t_iter,_) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
    }
    // Random-Scan Componentwise Estimation
    Rcpp::NumericVector u = log(runif(LIV));
    Rcpp::NumericVector z = rnorm(LIV);
    Rcpp::IntegerVector LIVseq = Rcpp::Range(1, LIV);
    Rcpp::IntegerVector s = Rcpp::RcppArmadillo::sample(LIVseq, LIV,
      false, NumericVector::create());
    s = s - 1;
    for (int j = 0; j < LIV; j++) {
      // Propose new parameter values
      Rcpp::List Mo0_ = clone(Mo0);
      Rcpp::NumericVector prop = Mo0_["parm"];
      prop[s[j]] += z[s[j]] * tuning[s[j]];
      // Log-Posterior of the proposed state
      Rcpp::List Mo1 = Model(prop, Data);
      fins = ::R_finite(Mo1["LP"]) + ::R_finite(Mo1["Dev"]);
      for (int m = 0; m < mcols; m++) {
        fins += ::R_finite(as<Rcpp::NumericVector>(Mo1["Monitor"])[m]);
      }
      if (fins < (mcols + 2)) Mo1 = Mo0;
      // Accept/Reject
      double LP0 = Mo0_["LP"];
      double LP1 = Mo1["LP"];
      log_alpha = LP1 - LP0;
      if (u[s[j]] < log_alpha) {
	Mo0 = Mo1;
        Acceptance[s[j]] += 1.0;
      }
    }
    if ((iter + 1) % Thinning == 0) {
      thinned(t_iter,_) = as<Rcpp::NumericVector>(Mo0["parm"]);
      Dev(t_iter,_) = as<Rcpp::NumericVector>(Mo0["Dev"]);
      Mon(t_iter,_) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
    }
    // Adapt the Proposal Variance
    if ((iter + 1) % Periodicity == 0) {
      double sqiter = pow(iter + 1.0, 0.5);
      if (sqiter > 100) sqiter = 100.0;
      double size = 1.0 / sqiter;
      Rcpp::NumericVector AcceptanceRate(Rcpp::clone(Acceptance));
      AcceptanceRate = AcceptanceRate / iter;
      Rcpp::NumericVector logtuning(log(tuning));
      for (int k = 0; k < LIV; k++) {
        if (AcceptanceRate[k] > 0.44) logtuning[k] += size;
        else logtuning[k] -= size;
	tuning = exp(logtuning);
        a_iter = floor(iter / Periodicity);
        DiagCovar(a_iter,_) = tuning;
      }
    }
  }
  int sumAccept = 0;
  for (int k = 0; k < LIV; k++) {
    sumAccept += Acceptance[k];
  }
  double acceptance = sumAccept / (LIV * 1.0);
  return wrap(Rcpp::List::create(Rcpp::Named("Acceptance") = acceptance,
                            Rcpp::Named("Dev") = Dev,
                            Rcpp::Named("DiagCovar") = DiagCovar,
                            Rcpp::Named("Mon") = Mon,
                            Rcpp::Named("thinned") = thinned,
                            Rcpp::Named("VarCov") = tuning));
}

SEXP mcmcmwg(SEXP MODEL, SEXP DATA, SEXP ITERATIONS, SEXP STATUS, 
  SEXP THINNING, SEXP SPECS, SEXP ACCEPTANCE, SEXP DEV, SEXP DIAGCOVAR, 
  SEXP liv, SEXP MON, SEXP MO0, SEXP SCALEF, SEXP THINNED, SEXP TUNING) {
  Function Model(MODEL);
  Rcpp::List Data(DATA), Specs(SPECS), Mo0(MO0);
  int Iterations = as<int>(ITERATIONS), Status = as<int>(STATUS),
    Thinning = as<int>(THINNING), LIV = as<int>(liv), 
    ScaleF = as<int>(SCALEF);
  double Acceptance = as<double>(ACCEPTANCE);
  Rcpp::NumericMatrix Dev(DEV), DiagCovar(DIAGCOVAR), Mon(MON),
    thinned(THINNED);
  Rcpp::NumericVector tuning(TUNING);
  int t_iter = 0, fins = 0, mcols = Mon.ncol();
  double log_alpha = 0;
  Rcpp::List Mo1 = clone(Mo0);
  RNGScope scope;
  for (int iter = 0; iter < Iterations; iter++) {
    // Print Status
    if ((iter + 1) % Status == 0) {
      Rcpp::Rcout << "Iteration: " << iter + 1 << 
           ",   Proposal: Componentwise,   LP: " << 
           floor(as<double>(Mo0["LP"]) * 100) / 100 << std::endl;
    }
    // Save Thinned Samples
    if ((iter + 1) % Thinning == 0) {
      t_iter = floor((iter) / Thinning) + 1;
      thinned(t_iter,_) = as<Rcpp::NumericVector>(Mo0["parm"]);
      Dev(t_iter,_) = as<Rcpp::NumericVector>(Mo0["Dev"]);
      Mon(t_iter,_) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
    }
    // Random-Scan Componentwise Estimation
    Rcpp::NumericVector u = log(runif(LIV));
    Rcpp::NumericVector z = rnorm(LIV);
    Rcpp::IntegerVector LIVseq = Rcpp::Range(1, LIV);
    Rcpp::IntegerVector s = Rcpp::RcppArmadillo::sample(LIVseq, LIV,
      false, NumericVector::create());
    s = s - 1;
    for (int j = 0; j < LIV; j++) {
      // Propose new parameter values
      Rcpp::List Mo0_ = clone(Mo0);
      Rcpp::NumericVector prop = Mo0_["parm"];
      prop[s[j]] += z[s[j]] * tuning[s[j]];
      // Log-Posterior of the proposed state
      Rcpp::List Mo1 = Model(prop, Data);
      fins = ::R_finite(Mo1["LP"]) + ::R_finite(Mo1["Dev"]);
      for (int m = 0; m < mcols; m++) {
        fins += ::R_finite(as<Rcpp::NumericVector>(Mo1["Monitor"])[m]);
      }
      if (fins < (mcols + 2)) Mo1 = Mo0;
      // Accept/Reject
      double LP0 = Mo0_["LP"];
      double LP1 = Mo1["LP"];
      log_alpha = LP1 - LP0;
      if (u[s[j]] < log_alpha) {
	Mo0 = Mo1;
        Acceptance += 1.0 / LIV;
      }
    }
    if ((iter + 1) % Thinning == 0) {
      thinned(t_iter,_) = as<Rcpp::NumericVector>(Mo0["parm"]);
      Dev(t_iter,_) = as<Rcpp::NumericVector>(Mo0["Dev"]);
      Mon(t_iter,_) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
    }
  }
  return wrap(Rcpp::List::create(Rcpp::Named("Acceptance") = Acceptance,
                            Rcpp::Named("Dev") = Dev,
                            Rcpp::Named("DiagCovar") = DiagCovar,
                            Rcpp::Named("Mon") = Mon,
                            Rcpp::Named("thinned") = thinned,
                            Rcpp::Named("VarCov") = tuning));
}
