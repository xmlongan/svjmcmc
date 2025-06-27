#include <Rcpp.h>
using namespace Rcpp;

double crsigma_v(double sigma_v_old, double prior_shape, double prior_rate,
                 NumericVector v, NumericVector y,
                 NumericVector jsize, NumericVector jtime,
                 double mu, double k,
                 double theta, double rho, double h);
