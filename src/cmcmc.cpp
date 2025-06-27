#include <Rcpp.h>
#include <cmath>

#include "crmu.h"
#include "crk.h"
#include "crtheta.h"
#include "crsigma_v.h"
#include "crrho.h"
#include "crlambda.h"
#include "crsigma_j.h"

#include "crv0.h"
#include "crvi.h"
#include "crvN.h"

#include "crjsize_i.h"
#include "crjtime_i.h"

#include "cinitialize_values.h"

using namespace Rcpp;

//' @rdname mcmc
// [[Rcpp::export]]
NumericVector cmcmc(NumericVector y, NumericVector ini_par,
                     int g = 5000, int G = 10000, double h = 1,
                     double echo = 0) {
   // g: warm-up samples
   // G: total samples
   // NumericVector parameters (5);
   //
   double mu=ini_par[0], k=ini_par[1], theta=ini_par[2];
   double sigma_v=ini_par[3], rho=ini_par[4];
   double lmbd = ini_par[5], sigma_j= ini_par[6];
   //
   double prior_shape_v0, prior_rate_v0;
   prior_shape_v0 = 2*k*theta/std::pow(sigma_v, 2);
   prior_rate_v0  = 2*k/std::pow(sigma_v, 2);
   //
   int N = y.size();
   NumericVector v (N+1), jsize (N), jtime (N);
   //
   List v_jsize_jtime = cinitialize_values(y,ini_par,1);
   v = v_jsize_jtime[0];
   jsize = v_jsize_jtime[1];
   jtime = v_jsize_jtime[2];
   //
   double cum_mu=0.0, cum_k=0.0, cum_theta=0.0, cum_sigma_v=0.0, cum_rho=0.0;
   double cum_lmdb=0.0, cum_sigma_j=0.0;
   //
   double prior_mu, prior_var, prior_shape, prior_rate;
   for (int i=0; i < G; i++) {
     // update parameters
     //
     prior_mu  = 0;
     prior_var = 0.1*0.1;
     mu = crmu(prior_mu, prior_var, v, y, jsize, jtime, k, theta, sigma_v, rho, h);
     //
     prior_mu  = 0.01;
     prior_var = 1*1;
     k = crk(prior_mu, prior_var, v, y, jsize, jtime, mu, theta, sigma_v, rho, h);
     //
     prior_mu  = 0.01;
     prior_var = 1*1;
     theta = crtheta(prior_mu, prior_var, v, y, jsize, jtime, mu, k, sigma_v, rho, h);
     //
     prior_shape = 2;
     prior_rate  = 1;
     sigma_v = crsigma_v(sigma_v, prior_shape, prior_rate,v,y,jsize,jtime,mu,k,theta,rho,h);
     //
     rho = crrho(rho, v, y, jsize, jtime, mu, k, theta, sigma_v, h);
     //
     double prior_a = 2.0, prior_b= 40.0;
     lmbd = crlambda(prior_a, prior_b, jtime);
     prior_shape= 5.0;
     double prior_scale = 20.0;
     sigma_j = crsigma_j(prior_shape, prior_scale, jsize);
     // update v0
     v[0] = crv0(v[0], prior_shape_v0, prior_rate_v0, v[1], y[0],
                 jsize[0], jtime[0],
                 mu, k, theta, sigma_v, rho, h);
     // update v1:v_{N-1}
     for (int n=1; n < N; n++) {
       v[n] = crvi(v[n], v[n-1], v[n+1], y[n-1], y[n],
                   jsize[n-1], jsize[n], jtime[n-1], jtime[n],
                   mu,k,theta,sigma_v,rho,h);
       //
       jsize[n-1] = crjsize_i(jtime[n-1], v[n-1], v[n], y[n-1], mu, k, theta,
                              sigma_v, rho, lmbd, sigma_j, h);
       jtime[n-1] = crjtime_i(jsize[n-1], v[n-1], v[n], y[n-1], mu, k, theta,
                              sigma_v, rho, lmbd, sigma_j, h);
     }
     // update vN
     v[N] = crvN(v[N], v[N-1], y[N-1], jsize[N-1], jtime[N-1],
                 mu, k, theta, sigma_v, rho, h);
     jsize[N-1] = crjsize_i(jtime[N-1], v[N-1], v[N], y[N-1], mu, k, theta,
                            sigma_v, rho, lmbd, sigma_j, h);
     jtime[N-1] = crjtime_i(jsize[N-1], v[N-1], v[N], y[N-1], mu, k, theta,
                            sigma_v, rho, lmbd, sigma_j, h);
     //
     // cumulate parameters
     if (i > g-1) {
       cum_mu += mu;
       cum_k  += k;
       cum_theta += theta;
       cum_sigma_v += sigma_v;
       cum_rho += rho;
       cum_lmdb += lmbd;
       cum_sigma_j += sigma_j;
     }
     //
     if (echo == 1 && (i < 10 || i > G-11)) {
       Rcout << i+1 << "th:\nparameters: " << mu << "," << k << ",";
       Rcout << theta << "," << sigma_v << "," << rho << "\n";
       Rcout << "first 10 of v: " << v[0] << "," << v[1] << "," << v[2] << ",";
       Rcout << v[3] << "," << v[4] << "," << v[5] << "," << v[6] << ",";
       Rcout << v[7] << "," << v[8] << "," << v[9] << "\n";
     }
   }
   // estimate
   int n = G - g;
   NumericVector estimate={cum_mu/n,cum_k/n,cum_theta/n,cum_sigma_v/n,cum_rho/n,
                           cum_lmdb/n, cum_sigma_j/n};
   return estimate;
 }
