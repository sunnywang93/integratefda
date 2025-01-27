% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgp.R
\name{bm_kl_rd}
\alias{bm_kl_rd}
\title{Generate brownian motion using the KL decomposition on a random grid of points}
\usage{
bm_kl_rd(k, x, lambda_rate, norm_constant, norm_factor, xi_dist, ...)
}
\arguments{
\item{k}{Numeric, number of basis functions.}

\item{x}{Vector containing the sampling points.}

\item{lambda_rate}{Numeric, the polynomial rate of decay of eigenvalues.}

\item{norm_constant}{Numeric, a normalising constant to make the scores centered
when using distribution that doesn't have zero mean.}

\item{norm_factor}{Numeric, a normalising factor to make the scores unit variance
when using a distribution that doesn't have unit variance by default.}

\item{xi_dist}{Function, indicating the distribution in which the scores
should be simulated from. Arguments associated to this function should be passed as
separate arguments.}
}
\value{
List, containing the sampling points, observed points, and normalised scores.
}
\description{
Generate brownian motion using the KL decomposition on a random grid of points
}
