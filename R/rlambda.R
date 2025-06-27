#' Update \eqn{\lambda}
#'
#' @description
#' Parameter \eqn{\lambda} has a Beta posterior if its prior is also a Beta.
#'
#' @param prior_a the first conjugate prior shape
#' @param prior_b the second conjugate prior shape
#' @param jtime vector of jump indicators \[\eqn{I_1, ..., I_N}\], I_n = 1 or 0.
#'
#' @return new value for \eqn{\lambda}, a scale value.
#' @export
#'
#' @examples
#' prior_a = 2
#' prior_b = 40
#' jtime = rbinom(4, size=1, prob=0.01)
#' rlambda(prior_a, prior_b, jtime)
rlambda <- function(prior_a, prior_b, jtime) {
  N = length(jtime)
  n = sum(jtime) # number of total jumps
  a = prior_a + n
  b = prior_b + N - n
  return(rbeta(1, a, b))
}
