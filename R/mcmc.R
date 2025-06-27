#' Estimate the parameters through MCMC
#'
#' @description
#' Given \eqn{[y_1,\cdots, y_N]}, to estimate parameters of Heston Stochastic
#' Volatility Model using MCMC.
#' * `cmcmc()` implemented in C++, returns only a vector of the estimated
#' parameters.
#' * `mcmc()` implemented in R.
#'
#' @details
#' Priors for and updating functions of these parameters:
#'
#' |prior|updating function|
#' |:---|:---|
#' |\eqn{\mu \sim \mathcal{N}(0,0.1^2) }|  [rmu()]|
#' |\eqn{ k \sim \mathcal{N}(0.01,1^2) }|  [rk()]|
#' |\eqn{\theta \sim \mathcal{N}(0.01,1^2)}|  [rtheta()]|
#' |\eqn{\sigma_v \sim Gamma(shape=2,rate=1)}|  [rsigma_v()]|
#' |\eqn{\rho \sim Uniform(-1,1)}|  [rrho()]|
#' |\eqn{\lambda \sim Beta(2,40)}|  [rlambda()]|
#' |\eqn{\sigma_j^2 \sim InverseGamma(5.0,20)}| [rsigma_j()]|
#'
#' Prior for the starting volatility: an uniform distribution over a large
#' range. (Depreciated
#' \eqn{v_0 \sim Gamma(shape=2k\theta/\sigma_v^2, rate=2k/\sigma_v^2)}.)
#'
#' Updating functions of volatility:
#' * \eqn{v_0}: [rv0()]
#' * \eqn{v_n}: [rvi()] for \eqn{1\le n < N}
#' * \eqn{v_N}: [rvN()]
#'
#' @param y vector of returns \[\eqn{y_1, ..., y_N}\].
#' @param ini_par vector of parameter initial values.
#' @param g number of warm-up iteration, defaults to 5,000.
#' @param G number of total iteration, defaults to 10,000.
#' @param h time unit, defaults to 1.
#' @param echo TRUE or FALSE, whether or not echo the parameters and first 10
#' volatilities in the first and last ten iterations, defaults to FALSE.
#'
#' @return a list of two elements:
#' * a vector of the estimated parameters, \eqn{(\mu,k,\theta,\sigma_v,\rho,\lambda,\sigma_j)},
#' * a matrix record of volatility in which each row represents an
#' iteration result, i.e., with dimension of \eqn{G\times(N+1)}.
#' @export
#'
#' @examples
#' y = rep(0.125,20)
#' ini_par=c(0.0625,0.05,0.125,0.05,-0.35, 0.01, 0.5)
#' mcmc(y, ini_par)
mcmc <- function(y, ini_par, g=5000, G=10000, h=1, echo=FALSE) {
  # g: warm-up samples
  # G: total samples
  init_values = initialize_values(y, ini_par=ini_par, h=1)
  #
  mu      = init_values$parameters[1]
  k       = init_values$parameters[2]
  theta   = init_values$parameters[3]
  sigma_v = init_values$parameters[4]
  rho     = init_values$parameters[5]
  lmbd    = init_values$parameters[6]
  sigma_j = init_values$parameters[7]
  #
  prior_shape_v0 = 2*k*theta/sigma_v^2
  prior_rate_v0 = 2*k/sigma_v^2
  #
  v = init_values$v
  jsize = init_values$jsize
  jtime = init_values$jtime
  #
  N = length(y)
  #
  record_mu = rep(0,G)
  record_k  = rep(0,G)
  record_theta = rep(0,G)
  record_sigma_v = rep(0,G)
  record_rho = rep(0,G)
  record_lmbd = rep(0,G)
  record_sigma_j = rep(0,G)
  record_v = matrix(rep(0,G*(N+1)), nrow=G, ncol=N+1)
  #
  for (i in 1:G) {
    # update parameters
    #
    prior_mu = 0; prior_var = 0.1^2
    mu = rmu(prior_mu, prior_var, v, y, jsize, jtime, k, theta, sigma_v, rho, h)
    #
    prior_mu = 0.01; prior_var = 1^2
    k = rk(prior_mu, prior_var, v, y, jsize, jtime, mu, theta, sigma_v, rho, h)
    #
    prior_mu = 0.01; prior_var = 1^2
    theta = rtheta(prior_mu, prior_var, v, y, jsize, jtime, mu, k, sigma_v, rho, h)
    #
    prior_shape = 2; prior_rate = 1
    sigma_v = rsigma_v(sigma_v, prior_shape, prior_rate, v,y,jsize,jtime,mu,k,theta,rho,h)
    #
    rho = rrho(rho, v, y, jsize, jtime, mu, k, theta, sigma_v, h)
    #
    prior_a = 2; prior_b = 40
    lmbd = rlambda(prior_a, prior_b, jtime)
    #
    prior_shape = 5.0; prior_scale = 20
    sigma_j = rsigma_j(prior_shape, prior_scale, jsize)
    # update v0
    v[1] = rv0(v[1], prior_shape_v0, prior_rate_v0, v[2], y[1],
               jsize[1], jtime[1],
               mu, k, theta, sigma_v, rho, h)
    record_v[i,1] = v[1]
    # update  v1:v_{N-1}, jsize_1:jsize_{N-1}
    for (n in 1:(N-1)) {
      v[n+1] = rvi(v[n+1], v[n], v[n+2], y[n], y[n+1],
                   jsize[n], jsize[n+1], jtime[n], jtime[n+1],
                   mu, k, theta, sigma_v, rho, h)
      record_v[i,n+1] = v[n+1]
      #
      jsize[n] = rjsize_i(jtime[n], v[n], v[n+1], y[n], mu, k, theta, sigma_v,
                          rho, lmbd, sigma_j, h)
      jtime[n] = rjtime_i(jsize[n], v[n], v[n+1], y[n], mu, k, theta, sigma_v,
                          rho, lmbd, sigma_j, h)
    }
    # update vN
    v[N+1] = rvN(v[N+1], v[N], y[N], jsize[N], jtime[N], mu, k, theta, sigma_v,
                 rho, h)
    record_v[i,N+1] = v[N+1]
    # update jsize
    jsize[N] = rjsize_i(jtime[N], v[N], v[N+1], y[N], mu, k, theta, sigma_v,
                        rho, lmbd, sigma_j, h)
    jtime[N] = rjtime_i(jsize[N], v[N], v[N+1], y[N], mu, k, theta, sigma_v,
                        rho, lmbd, sigma_j, h)
    # record
    record_mu[i] = mu
    record_k[i] = k
    record_theta[i] = theta
    record_sigma_v[i] = sigma_v
    record_rho[i] = rho
    record_lmbd[i] = lmbd
    record_sigma_j[i] = sigma_j
    #
    if (echo==TRUE && (i<11 || i>G-10)) {
      cat(sprintf('\n%dth: \nparameters: ',i))
      cat(signif(c(mu,k,theta,sigma_v,rho), digits=4), sep=',')
      cat('\nfrist ten of v: ')
      cat(signif(v[1:10], digits = 3), sep = ',')
    }
  }
  #
  # estimate
  mu = mean(record_mu[g:G])
  k  = mean(record_k[g:G])
  theta = mean(record_theta[g:G])
  sigma_v = mean(record_sigma_v[g:G])
  rho = mean(record_rho[g:G])
  lmbd = mean(record_lmbd[g:G])
  sigma_j = mean(record_sigma_j[g:G])
  #
  return(list(parameters = c(mu,k,theta,sigma_v,rho,lmbd,sigma_j),
              record_v = record_v[seq(G-49,G), seq(N-9,N)]))
}
