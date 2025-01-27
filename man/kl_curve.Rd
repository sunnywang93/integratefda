% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/curves.R
\name{kl_curve}
\alias{kl_curve}
\title{Recover the curves using the KL-decomposition}
\usage{
kl_curve(scores_vec, psi_list)
}
\arguments{
\item{scores_vec}{Vector of scores}

\item{psi_list}{List of the eigenfunctions, containing the elements:
\itemize{
\item \strong{$t} Vector of sampling points.
\item \strong{$X} Matrix of observed points, with sampling points on rows and the index
of the eigenfunctions on the columns.
}}
}
\value{
List containing the sampling points and observed points.
}
\description{
Recover the curves using the KL-decomposition
}
