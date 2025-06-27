#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname rsigma_j
// [[Rcpp::export]]
double crsigma_j(double prior_shape, double prior_scale, NumericVector jsize) {
  int N = jsize.size();
  double shape = prior_shape + N/2.0;
  double scale = prior_scale + Rcpp::sum(jsize * jsize) / 2.0;
  Environment pkg = Environment::namespace_env("invgamma");
  Function rinvgamma = pkg["rinvgamma"];
  NumericVector sigma_j2_vec = rinvgamma(1, Named("shape", shape), Named("scale", scale));
  double sigma_j2 = sigma_j2_vec[0];
  return std::sqrt(sigma_j2);
}
