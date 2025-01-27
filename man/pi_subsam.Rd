% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mc_nn.R
\name{pi_subsam}
\alias{pi_subsam}
\title{Estimates prediction intervals using subsampling}
\usage{
pi_subsam(varphi, varphi_int, eps = 0.05, s, b_out = 200, int_fun)
}
\arguments{
\item{varphi}{Vector, the integrand on a vector of evaluation points.}

\item{varphi_int}{Numeric, the integral estimated.}

\item{eps}{Numeric, the error critical value.}

\item{s}{Numeric, the Hölder exponent of the integrand.}

\item{b_out}{Numeric, number of bootstrap replications to perform.}

\item{int_fun}{Function, indicating the method of approximating integrals.
Only \code{mean} and \code{trapz} function is accepted.}
}
\value{
List, containing the following elements:
\itemize{
\item \strong{$varphi_int} Numeric, the estimated integral value.
\item \strong{$pi_l} Numeric, the lower bound of the prediction interval.
\item \strong{$pi_u} Numeric, the upper bound of the prediction interval.
\item \strong{$width} Numeric, the width of the prediction interval.
}
}
\description{
Mainly used for comparison purposes. Only works with either the sample mean
or riemann sums using the trapezoidal rule.
}
