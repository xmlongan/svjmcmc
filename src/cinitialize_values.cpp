#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname initialize_values
// [[Rcpp::export]]
List cinitialize_values(NumericVector y, NumericVector ini_par, double h) {
   double std_eps_y, mu_v, sd_v;
   // initialize parameter values
   // mu=0, k=0.01, theta=0.1, sigma_v=0.01, rho=0
   // parameters[0] = 0;
   // parameters[1] = 0.01;
   // parameters[2] = 0.1;
   // parameters[3] = 0.01;
   // parameters[4] = 0;
   // double mu=0, k=0.01, theta=0.1, sigma_v=0.01, rho=0;
   double mu=ini_par[0], k=ini_par[1], theta=ini_par[2];
   double sigma_v=ini_par[3], rho=ini_par[4];
   double lmbd = ini_par[5], sigma_j= ini_par[6];
   //
   int N = y.size();
   NumericVector v (N+1);
   // initial jsize, jtime
   NumericVector jsize (N), jtime (N);
   jsize = rnorm(N, 0.0, sigma_j);
   jtime = rbinom(N, 1.0, lmbd*h);
   // let v_0 be the long-run mean of volatility
   v[0] = theta;
   //
   for (int n=1; n < N+1; n++) {
     // y_n <- v_{n-1}
     std_eps_y = y[n-1] - mu*h + v[n-1]*h/2 - jsize[n-1]*jtime[n-1];
     // v_n <- v_{n-1},y_n
     mu_v = v[n-1] + k * (theta - v[n-1]) * h + sigma_v * rho * std_eps_y;
     sd_v = sigma_v * std::sqrt((1-rho*rho) * v[n-1] * h);
     v[n] = R::rnorm(mu_v, sd_v);
     if (v[n] <= 0) {
       v[n] = 0.00001;
     }
   }
   return List::create(v, jsize, jtime);
 }
