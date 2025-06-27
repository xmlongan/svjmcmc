#' Update \eqn{\sigma_y^2}
#'
#' @description
#' Parameter \eqn{\sigma_y^2} has an Inverse Gamma posterior if its prior is
#' also an Inverse Gamma.
#'
#' @param prior_shape conjugate prior shape parameter
#' @param prior_scale conjugate prior scale parameter
#' @param jsize vector of jump sizes \[\eqn{J_1, ..., J_N}\]
#'
#' @return new value for \eqn{\sigma_y}, a scale value.
#' @export
#'
#' @examples
#' prior_shape = 5.0
#' prior_scale = 20
#' jsize = rnorm(4, mean=0, sd=0.5)
#' rsigma_j(prior_shape, prior_scale, jsize)
rsigma_j <- function(prior_shape, prior_scale, jsize) {
  N = length(jsize)
  shape = prior_shape + N/2
  scale = prior_scale + sum(jsize^2)/2. # mu_j = 0
  sigma_j2 = invgamma::rinvgamma(1, shape=shape, scale=scale)
  sigma_j = sqrt(sigma_j2)
  return(sigma_j)
}
