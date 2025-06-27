#' Update \eqn{J_i}
#'
#' @description
#' Generate new jump size.
#'
#' @param jtime_i jump indicator at epoch i.
#' @param vim1 variance \eqn{v_{i-1}}.
#' @param vi variance \eqn{v_i}.
#' @param yi return \eqn{y_i}.
#' @param mu parameter \eqn{\mu}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma_v parameter \eqn{\sigma_v}.
#' @param rho parameter \eqn{\rho}.
#' @param lmbd parameter \eqn{\lambda}.
#' @param sigma_j paramter \eqn{\sigma_j}.
#' @param h time unit.
#'
#' @return new value of \eqn{J_i}, a scale value.
#' @export
#'
rjsize_i <- function(jtime_i, vim1, vi, yi, mu,k,theta,sigma_v,rho,lmbd,sigma_j,h) {
  if (jtime_i == 0) {
    return(rnorm(1, mean=0, sd=sigma_j))
  }
  # now jtime_i == 1
  std_eps_v = (vi - vim1 - k*(theta-vim1)*h)/sigma_v
  den = (1-rho^2)*vim1*h + sigma_j^2

  mean = sigma_j^2 * (yi - mu*h + vim1*h/2 - rho*std_eps_v) / den
  sd = sigma_j * sqrt((1-rho^2)*vim1*h / den)
  return(rnorm(1, mean=mean, sd=sd))
}
