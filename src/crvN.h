#include <Rcpp.h>
using namespace Rcpp;

double crvN(double vN_old, double vNm1, double yN, double jsize_N, double jtime_N,
            double mu, double k,
            double theta, double sigma_v, double rho, double h);
