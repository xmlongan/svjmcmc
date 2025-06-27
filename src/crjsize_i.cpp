#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname rjsize_i
// [[Rcpp::export]]
double crjsize_i(double jtime_i, double vim1, double vi, double yi,
                 double mu, double k, double theta, double sigma_v, double rho,
                 double lmbd, double sigma_j, double h) {
  if (jtime_i == 0.0) {
    return R::rnorm(0.0, sigma_j);
  }
  // now jtime_i = 1
  double std_eps_v = (vi - vim1 - k*(theta-vim1)*h)/sigma_v;
  double den = (1-rho*rho)*vim1*h + sigma_j*sigma_j;

  double mean = sigma_j*sigma_j * (yi - mu*h + vim1*h/2 - rho*std_eps_v)/den;
  double sd = sigma_j * std::sqrt((1-rho*rho)*vim1*h/den);
  return R::rnorm(mean, sd);
}
