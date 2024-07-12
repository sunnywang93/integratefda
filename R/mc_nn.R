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

  list(
    t = wm_loo$x,
    weights = wm_loo$weights,
    varphi_int = sum(wm_loo$weights * varphi)
  )

}

#' Computes confidence intervals with control neighbours integration method in
#' the presence of measurement errors
#'
#' @param varphi_noise List, containing the following elements:
#' - **$weights** Vector of control neighbour weights.
#' - **$varphi_int** Numeric, the estimated integral value.
#' - **$varphi_noise** Vector, the noise functional in which the
#' confidence intervals should be based on.
#' @param conf_level Numeric, the desired confidence level.
#' @param lim Boolean, indicating whether to use the asymptotic confidence intervals.
#' If FALSE, the non-asymptotic version based on conditional variance is used.
#' @returns List, containing the following elements:
#' - **$varphi_int** Numeric, the estimated integral value.
#' - **$ci_l** Numeric, the lower bound of the confidence interval.
#' - **$ci_u** Numeric, the upper bound of the confidence interval.
#' @export
#'
mc_ci <- function(varphi_noise, conf_level = 0.95, lim = FALSE) {

  z_value <- qnorm(1 - (1 - conf_level)/2)

  if(lim) {
    sm <- sqrt(5/2) * length(varphi_noise$weights)^(-1/2) *
      mean(varphi_noise$varphi_noise^2)^(1/2)
    ci_l <- varphi_noise$varphi_int - z_value * sm
    ci_u <- varphi_noise$varphi_int + z_value * sm

  } else {
    sm <- sqrt(sum(varphi_noise$weights^2 * varphi_noise$varphi_noise^2))
    ci_l <- varphi_noise$varphi_int - z_value * sm
    ci_u <- varphi_noise$varphi_int + z_value * sm
  }


  list(varphi_int = varphi_noise$varphi_int,
       ci_l = ci_l,
       ci_u = ci_u,
       sm = sm,
       width = abs(ci_u - ci_l))

}






