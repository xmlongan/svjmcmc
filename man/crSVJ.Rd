% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{crSVJ}
\alias{crSVJ}
\title{Generate SVJ Sample Trajectory}
\usage{
crSVJ(v_0, n_segment, par, N, h = 1)
}
\arguments{
\item{v_0}{the initial volatility \eqn{v_0}.}

\item{n_segment}{number of further segments for each interval.}

\item{par}{a vector of true values for the parameters
\eqn{\mu,k,\theta,\sigma_v,\rho,\lambda,\mu_j,\sigma_j}.}

\item{N}{the number of sample returns, \eqn{y_1,\cdots,y_N}.}

\item{h}{time unit.}
}
\value{
a vector of returns, \eqn{(y_1,\cdots,y_N)}.
}
\description{
Generate SVJ sample trajectory using Euler
Approximation.
}
\details{
SVJ:
\deqn{y_n=\mu h -\frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}(\rho\epsilon_n^v
      +\sqrt{1-\rho^2}\epsilon_n) + J_nI_n,}
\deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}\epsilon_n^v.}
}
\examples{
S0 = c(0.125,0.1,0.25,0.1,-0.7, 0.1, 0.0, 0.5)
y_series_0 = crSVJ(v_0=S0[3], n_segment=10, par=S0, N=1000, h=1)

S1 = c(0.4,0.1,0.25,0.1,-0.7, 0.1, 0.0, 0.5)
y_series_1 = crSVJ(v_0=S1[3], n_segment=10, par=S1, N=1000, h=1)

S2 = c(0.125,0.03,0.25,0.1,-0.7, 0.1, 0.0, 0.5)
y_series_2 = crSVJ(v_0=S2[3], n_segment=10, par=S2, N=1000, h=1)

S3 = c(0.125,0.1,0.5,0.1,-0.7, 0.1, 0.0, 0.5)
y_series_3 = crSVJ(v_0=S3[3], n_segment=10, par=S3, N=1000, h=1)

S4 = c(0.125,0.1,0.25,0.2,-0.7, 0.1, 0.0, 0.5)
y_series_4 = crSVJ(v_0=S4[3], n_segment=10, par=S4, N=1000, h=1)

S5 = c(0.125,0.1,0.25,0.1,-0.3, 0.1, 0.0, 0.5)
y_series_5 = crSVJ(v_0=S5[3], n_segment=10, par=S5, N=1000, h=1)
}
