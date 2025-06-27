#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname rjtime_i
// [[Rcpp::export]]
double crjtime_i(double jsize_i, double vim1, double vi, double yi, double mu,
                  double k, double theta, double sigma_v, double rho,
                  double lmbd, double sigma_j, double h) {
   double std_eps_v = (vi - vim1 - k*(theta-vim1)*h)/sigma_v;
   double xi = yi - mu*h + vim1*h/2.0 -rho*std_eps_v;
   double den = 2.0*(1-rho*rho)*vim1*h;
   double weight = std::exp((jsize_i*jsize_i - 2*xi*jsize_i)/den);

   double post_prob = lmbd*h / (lmbd*h + (1-lmbd*h)*weight);
   return R::rbinom(1.0, post_prob);
 }
