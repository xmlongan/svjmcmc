#include <Rcpp.h>
using namespace Rcpp;

double crvi(double vi_old, double vim1, double vip1, double yi, double yip1,
            double jsize_i, double jsize_ip1, double jtime_i, double jtime_ip1,
            double mu, double k, double theta, double sigma_v, double rho,
            double h);
