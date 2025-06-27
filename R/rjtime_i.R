#' Update \eqn{I_i}
#'
#' @description
#' Generate new jump indicator, whether there is a jump or not.
#'
#' @param jsize_i jump size at epoch i.
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
rjtime_i <- function(jsize_i, vim1, vi, yi, mu,k,theta,sigma_v,rho,lmbd,sigma_j,h) {
  std_eps_v = (vi - vim1 - k*(theta-vim1)*h)/sigma_v
  xi <- yi - mu*h + vim1*h/2 - rho*std_eps_v
  den <- 2*(1-rho^2)*vim1*h
  weight <- exp((jsize_i^2 - 2*xi*jsize_i)/den)

  post_prob = lmbd*h / (lmbd*h + (1-lmbd*h) * weight)
  return(rbinom(1, size=1, prob=post_prob))
}
