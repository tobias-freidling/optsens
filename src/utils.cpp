#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat resid_cpp(const arma::mat& X, const arma::mat& Y) {
  return Y - X * arma::solve(X, Y);
}


// [[Rcpp::export]]
double r_marginal(const arma::vec& A, const arma::vec& B) {
  return arma::as_scalar(A.t() * B) / arma::norm(A,2) / arma::norm(B,2);
}


// [[Rcpp::export]]
double r_partial(const arma::vec& A, const arma::vec& B, const arma::mat& C) {
  
  arma::vec A_C = A - C * arma::solve(C, A);
  arma::vec B_C = B - C * arma::solve(C, B);
  
  return r_marginal(A_C, B_C);
}


// [[Rcpp::export]]
double r2_marginal(const arma::vec& A, const arma::mat& B) {
  arma::vec A_B = A - B * arma::solve(B, A);
  return 1 - arma::as_scalar(A_B.t() * A_B) / arma::as_scalar(A.t() * A);
}


// [[Rcpp::export]]
double r2_partial(const arma::vec& A, const arma::mat& B, const arma::mat& C) {
  
  arma::vec A_C = A - C * arma::solve(C, A);
  arma::mat B_C = B - C * arma::solve(C, B);
  
  return r2_marginal(A_C, B_C);
}
