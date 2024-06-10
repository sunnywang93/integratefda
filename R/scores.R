#' Compute scores using leave-one-out control neighbours monte carlo integration
#'
#' @param y_list List of the observations, containing the elements:
#' - **$t** Vector of sampling points.
#' - **$y** Vector of observed points.
#' @param mu_list List of the mean function, containing the elements:
#' - **$t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @param psi_list List of the eigenfunctions, containing the elements:
#' - **$t** Vector of sampling points.
#' - **$X** Matrix of observed points, with sampling points on rows and the index
#' of the eigenfunctions on the columns.
#' @param pdf_list List of the density, containing the elements:
#' - **$t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @param cdf List or a function representing the cumulative distribution function.
#' If list, it should contain the elements:
#' - **$t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @returns Vector containing the estimated scores.
#' @export


mc_scores <- function(y_list, mu_list, psi_list, pdf_list, cdf, noise_sd) {

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

  varphi <- (y_list$y - mu_list$x) * psi_list$X / pdf_list$x

  varphi_ci <- noise_sd * psi_list$X / pdf_list$x

  list(t = wm_loo$x,
       weights = wm_loo$weights,
       xi_hat = colSums(wm_loo$weights * varphi),
       varphi = varphi,
       varphi_ci = varphi_ci
  )

}


#' Computes the confidence intervals for scores based on the conditional variance
#' of leave-one-out control neighbours
#'
#' @param scores_list List, containing the elements:
#' - **$weights** Vector, containing the weights associated to each score.
#' - **$xi_hat** Vector, containing the estimated scores.
#' - **$varphi_ci** Matrix, containing the evaluated functional, with rows carrying
#' the sampling points and columns carrying the components.
#' @param conf_level Numeric, indicating the confidence level in which the
#' intervals should be constructed.
#' @returns List, containing the elements:
#' - **$scores** Vector, containing the estimated scores.
#' - **$ci_l** Vector, containing the lower bound of the confidence interval.
#' - **$ci_u** Vector, containing the upper bound of the confidence interval.
#' @export

confint_mc <- function(scores_list, conf_level) {

  sm <- sqrt(colSums(scores_list$weights^2 * scores_list$varphi_ci^2))

  crit <- 1 - conf_level

  ci_l <- scores_list$xi_hat - qnorm(1 - crit/2) * sm
  ci_u <- scores_list$xi_hat + qnorm(1 - crit/2) * sm

  list(scores = scores_list$xi_hat,
       ci_l = ci_l,
       ci_u = ci_u,
       sm = sm)

}

#' Computes the confidence intervals for scores based on the asymptotic limit
#' of leave-one-out control neighbours
#'
#' @param scores_list List, containing the elements:
#' - **$weights** Vector, containing the weights associated to each score.
#' - **$xi_hat** Vector, containing the estimated scores.
#' - **$varphi** Matrix, containing the evaluated functional, with rows carrying
#' the sampling points and columns carrying the components.
#' @param pdf_list List, containing the elements:
#' - **t** Vector containing the sampling points.
#' - **x** Vector containing the observed points.
#' @param conf_level Numeric, indicating the confidence level in which the
#' intervals should be constructed.
#' @returns List, containing the elements:
#' - **$scores** Vector, containing the estimated scores.
#' - **$ci_l** Vector, containing the lower bound of the confidence interval.
#' - **$ci_u** Vector, containing the upper bound of the confidence interval.
#' - **$sm** Vector, containing the standard deviation.
#' @export


confint_lim <- function(scores_list, pdf_list, conf_level) {

  int_varphi <- apply(scores_list$varphi_ci^2 * pdf_list$x,
                      2,
                      function(x) pracma::trapz(x = pdf_list$t,
                                                y = x)
  )


  sm <- sqrt(5/2) * (length(scores_list$weights))^(-1/2) * (int_varphi)^(1/2)

  crit <- 1 - conf_level

  ci_l <- scores_list$xi_hat - qnorm(1 - crit/2) * sm
  ci_u <- scores_list$xi_hat + qnorm(1 - crit/2) * sm

  list(scores = scores_list$xi_hat,
       ci_l = ci_l,
       ci_u = ci_u,
       sm = sm)


}




