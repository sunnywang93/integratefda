% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.R
\name{mc_scores}
\alias{mc_scores}
\title{Compute scores using leave-one-out control neighbours integration method.}
\usage{
mc_scores(y_list, mu_list, psi_list, pdf_list, cdf, noise_sd)
}
\arguments{
\item{y_list}{List of the observations, containing the elements:
\itemize{
\item \strong{$t} Vector of sampling points.
\item \strong{$y} Vector of observed points.
}}

\item{mu_list}{List of the mean function, containing the elements:
\itemize{
\item \strong{$t} Vector of sampling points.
\item \strong{$x} Vector of observed points.
}}

\item{psi_list}{List of the eigenfunctions, containing the elements:
\itemize{
\item \strong{$t} Vector of sampling points.
\item \strong{$X} Matrix of observed points, with sampling points on rows and the index
of the eigenfunctions on the columns.
}}

\item{pdf_list}{List of the density, containing the elements:
\itemize{
\item \strong{$t} Vector of sampling points.
\item \strong{$x} Vector of observed points.
}}

\item{cdf}{List or a function representing the cumulative distribution function.
If list, it should contain the elements:
\itemize{
\item \strong{$t} Vector of sampling points.
\item \strong{$x} Vector of observed points.
}}

\item{noise_sd}{Numeric or vector, representing the standard deviation of the
noise that can depend on the design points.}
}
\value{
Vector containing the estimated scores.
}
\description{
Compute scores using leave-one-out control neighbours integration method.
}
