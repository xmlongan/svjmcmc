% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R, R/rv0.R
\name{crv0}
\alias{crv0}
\alias{rv0}
\title{Update \eqn{v_0}}
\usage{
crv0(
  v0_old,
  prior_shape,
  prior_rate,
  v1,
  y1,
  jsize_1,
  jtime_1,
  mu,
  k,
  theta,
  sigma_v,
  rho,
  h
)

rv0(
  v0_old,
  prior_shape,
  prior_rate,
  v1,
  y1,
  jsize_1,
  jtime_1,
  mu,
  k,
  theta,
  sigma_v,
  rho,
  h
)
}
\arguments{
\item{v0_old}{value of \eqn{v_0} before updating, i.e.,
\eqn{v_0^{(g)}}, in the \eqn{g}th iteration.}

\item{prior_shape}{parameter shape of \eqn{v_0}'s prior gamma distribution
\eqn{P_{prior}(\cdot)}.}

\item{prior_rate}{parameter rate of \eqn{v_0}'s prior gamma distribution
\eqn{P_{prior}(\cdot)}.}

\item{v1}{value of \eqn{v_1}.}

\item{y1}{value of \eqn{y_1}.}

\item{jsize_1}{jump size \eqn{J_1}.}

\item{jtime_1}{jump or not \eqn{I_1}, 0 or 1.}

\item{mu}{parameter \eqn{\mu}.}

\item{k}{parameter \eqn{k}.}

\item{theta}{parameter \eqn{\theta}.}

\item{sigma_v}{parameter \eqn{\sigma_v}.}

\item{rho}{parameter \eqn{\rho}.}

\item{h}{time unit.}
}
\value{
new \eqn{v_0}, a scale value.
}
\description{
Generate new first volatility \eqn{v_0} according to its posterior,
using Independence
Metropolis with a gamma density proposal \eqn{P(\cdot)}, i.e.,
`rgamma(1, shape=2, rate=1/v0_old)`, whoes mode equals *v0_old*.
The prior of \eqn{v_0} is assumed to be an uniform over a very large range.
(*Depreciated*: Its prior is
a gamma distribution with shape \eqn{2k\theta/\sigma_v^2} and rate
\eqn{2k/\sigma_v^2}, *requiring* \eqn{k,\theta,\sigma_v} *at least being
positive*. Note that the shape and rate are computed and remain
unchanged, using the initial values of the parameters.)
* `crv0()` implemented in C++, through package Rcpp.
* `rv0()` implemented in R.
}
\section{Major steps}{


1. First, propose a new \eqn{v_0^{(g+1)}} according to gamma distribution,
`rgamma(1, shape=2, rate=1/v0_old)`.
2. Then, compute the proposal density ratio
   \deqn{\frac{P(v_0^{(g)})}{P(v_0^{(g+1)})}=
         \frac{dgamma(v_0^{(g)},shape=2,rate=1/v_0^{(g+1)})}
              {dgamma(v_0^{(g+1)},shape=2,rate=1/v_0^{(g)})}.}
3. Last, accept \eqn{v_0^{(g+1)}} with probability
\deqn{\alpha(v_0^{(g)},v_0^{(g+1)})
=min\left(\frac{\pi(v_0^{(g+1)})}{\pi(v_0^{(g)})}
         \times \frac{P(v_0^{(g)})}{P(v_0^{(g+1)})},1\right).}
}

\section{Posterior density ratio}{


Employ the following version of SVJ model
\deqn{y_n=\mu h -\frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}\epsilon_n^y + J_nI_n}
\deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}(
     \rho\epsilon_n^y + \sqrt{1-\rho^2}\epsilon_n)}
The posterior
\deqn{\begin{matrix}
       P(v_0|v_{1:N},y_{1:N},J_{1:N},I_{1:N})
       &\propto& P(v_0,v_{1:N},y_{1:N}|J_{1:N},I_{1:N})\\
       &\propto& P(v_1|v_0,y_1,J_1,I_1)\cdot P(y_1|v_0,J_1,I_1) \cdot P_{prior}(v_0)
      \end{matrix}}
Introduce notation
\deqn{\tilde{v}_1\triangleq v_0+k(\theta-v_0)h+\sigma_v\rho
      \sqrt{v_0h}\epsilon_1^y,\quad
      \sqrt{v_0h}\epsilon_1^y = y_1-\mu h +\frac{1}{2}v_0h - J_1I_1,}
thus, we have
\deqn{v_1|v_0,y_1,J_1,I_1
      \sim \mathcal{N}(\tilde{v}_1, \sigma_v^2(1-\rho^2)v_0h),\quad
      P(v_1|v_0,y_1,J_1,I_1)
      = \frac{1}{\sqrt{2\pi(1-\rho^2)v_0h}\sigma_v}
        e^{-\frac{(v_1-\tilde{v}_1)^2}{2\sigma_v^2(1-\rho^2)v_0h}}.}
The posterior density ratio
\deqn{\begin{matrix}
  \frac{\pi(v_0^{(g+1)})}{\pi(v_0^{(g)})}
    &=&\frac{P(v_1|v_0^{(g+1)},y_1)}{P(v_1|v_0^{(g)},y_1)}\times
       \frac{P(y_1|v_0^{(g+1)})}{P(y_1|v_0^{(g)})}\times
       \frac{P_{prior}(v_0^{(g+1)})}{P_{prior}(v_0^{(g)})}\\
  log\left(\frac{\pi(v_0^{(g+1)})}{\pi(v_0^{(g)})}\right)
    &=&-log\frac{v_0^{(g+1)}}{v_0^{(g)}}
   -\frac{(v_1-\mu^{(g+1)})^2}{2\sigma_v^2(1-\rho^2)hv_0^{(g+1)}}
   +\frac{(v_1-\mu^{(g)})^2}{2\sigma_v^2(1-\rho^2)hv_0^{(g)}}\\
   &&-\frac{(y_1-\mu h +v_0^{(g+1)}h/2)^2}{2hv_0^{(g+1)}}
   +\frac{(y_1-\mu h +v_0^{(g)}h/2)^2}{2hv_0^{(g)}}
   +log\frac{P_{prior}(v_0^{(g+1)})}{P_{prior}(v_0^{(g)})}
\end{matrix}}
where
\deqn{\mu^{(g)}\triangleq v_0^{(g)} + k(\theta-v_0^{(g)})h + \sigma_v\rho
\sqrt{v_0^{(g)}h}\epsilon_1^y}
substitute \eqn{g} with \eqn{g+1} to get \eqn{\mu^{(g+1)}}.
}

\examples{
v1 = 0.2; y1 = 0
mu = 0.01; k = 0.1; theta = 0.25; sigma_v = 0.25; rho = 0; h = 1
prior_shape = 2*k*theta/sigma_v^2
prior_rate = 2*k/sigma_v^2
jsize_1 = rnorm(1, mean=0, sd=0.5)
jtime_1 = rbinom(1,size=1, prob=0.01)
rv0(0.1,prior_shape,prior_rate,v1,y1,jsize_1,jtime_1,mu,k,theta,sigma_v,rho,h)
}
