#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname rlambda
// [[Rcpp::export]]
double crlambda(double prior_a, double prior_b, NumericVector jtime) {
  int N = jtime.size();
  int n = Rcpp::sum(jtime);
  double a = prior_a + n;
  double b = prior_b + N - n;
  return R::rbeta(a, b);
}
