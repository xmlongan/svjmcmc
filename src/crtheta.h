#include <Rcpp.h>
using namespace Rcpp;

double crtheta(double prior_mu, double prior_var, NumericVector v,
               NumericVector y, NumericVector jsize, NumericVector jtime,
               double mu, double k, double sigma_v, double rho,
               double h);
