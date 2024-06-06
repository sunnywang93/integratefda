#' Compute scores using leave-one-out control neighbours monte carlo integration
#'
#' @param y_list List of the observations, containing the elements:
#' -**$t** Vector of sampling points.
#' -**$y** Vector of observed points.
#' @param mu_list List of the mean function, containing the elements:
#' -**$t** Vector of sampling points.
#' -**$x** Vector of observed points.
#' @param psi_list List of the eigenfunctions, containing the elements:
#' -**$t** Vector of sampling points.
#' -**$X** Matrix of observed points, with sampling points on rows and the index
#' of the eigenfunctions on the columns.
#' @param pdf_list List of the density, containing the elements:
#' -**$t** Vector of sampling points.
#' -**$x** Vector of observed points.
#' @param cdf List or a function representing the cumulative distribution function.
#' If list, it should contain the elements:
#' -**$t** Vector of sampling points.
#' -**$x** Vector of observed points.
#' @returns Vector containing the estimated scores.
#' @export


mc_scores <- function(y_list, mu_list, psi_list, pdf_list, cdf) {

  if(!is.function(cdf)) {
    cdf_fun <- function(x) {
      id <- sapply(x, function(xi) which.min(abs(xi - cdf$t)))
      cdf$x[id]
    }
  } else {
    cdf_fun <- cdf
  }

  wm_loo <- weights_loo(x = y_list$t,
                        cdf = cdf_fun)

  colSums((wm_loo$weights * (y_list$y - mu_list$x) * psi_list$X) / pdf_list$x)

}





