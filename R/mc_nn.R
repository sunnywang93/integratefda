#' Estimates the integral of a function using control neighbours method
#'
#' @param x Vector of sampling points.
#' @param varphi Vector containing the observed points of the integrand.
#' @param cdf List or a function representing the cumulative distribution function.
#' If list, it should contain the elements:
#' - **$t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @returns Numeric, the estimated integral.
#' @export
#'


mc_int <- function(x, varphi, cdf) {

  if(!is.function(cdf)) {
    cdf_fun <- function(x) {
      id <- sapply(x, function(xi) which.min(abs(xi - cdf$t)))
      cdf$x[id]
    }
  } else {
    cdf_fun <- cdf
  }

  wm_loo <- weights_loo(x = x,
                        cdf = cdf_fun)

  sum(wm_loo$weights * varphi)

}









