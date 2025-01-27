% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mc_nn.R
\name{mc_pi2d}
\alias{mc_pi2d}
\title{Computes prediction intervals for 2D control neighbour integral estimates
based on subsampling}
\usage{
mc_pi2d(X, X_int, eps = 0.05, s, b_out = 200)
}
\arguments{
\item{X}{Data frame, containing the following columns:
\itemize{
\item \strong{$t1} Vector, containing the coordinates in \code{t1}.
\item \strong{$t2} Vector, containing the coordinates in \code{t2}.
\item \strong{$x} Vector of observed points, where each point is associated with each
row of coordinates in \code{t}.
}}

\item{X_int}{Numeric, the control neighbour integral estimate.}

\item{eps}{Numeric, the error critical value.}

\item{s}{Numeric, the Hölder exponent of the integrand.}

\item{b_out}{Numeric, number of bootstrap replications to perform.}
}
\value{
List, containing the following elements:
\itemize{
\item \strong{$X_int} Numeric, the estimated integral value.
\item \strong{$pi_l} Numeric, the lower bound of the prediction interval.
\item \strong{$pi_u} Numeric, the upper bound of the prediction interval.
\item \strong{$width} Numeric, the width of the prediction interval.
}
}
\description{
Computes prediction intervals for 2D control neighbour integral estimates
based on subsampling
}
