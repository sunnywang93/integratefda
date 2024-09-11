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

mc_pi_bound <- function(varphi, varphi_int, M, L, eps, s, b = 1,
                  cvol = 1/2, U = 1, d = 1) {

  if(eps <= 2 / M^(1/2 + s/d)) {
    delta_eps <- log(24 * M / eps)^(1 + s/d)
  } else {
    delta_eps <- (3*log(3*M))^(1 + s/d)
  }

  Delta <- (sqrt(log(2 / sqrt(eps))) + M^(-s/d)) / (M^(1/2 + s/d)) * delta_eps
  Cphi <- 1/2 * (max(varphi) - min(varphi))
  Cphi_norm <- Cphi / M^(1/2 + s/d)
  ckiss <- 2^(0.401*d)

  Vd <- pi^(d/2) / gamma(d/2 + 1)

  K1 <- 4 * L * ((17 * d) / (b * cvol * Vd))^(s/d) *
    (ckiss + 1 + ((17 * d * U) / (b*cvol)))

  ci_u <- varphi_int + (K1 * Delta + Cphi_norm)
  ci_l <- varphi_int - (K1 * Delta + Cphi_norm)

  list(varphi_int = varphi_int,
       ci_u = ci_u,
       ci_l = ci_l,
       width = abs(ci_u - ci_l))

}


#' Computes prediction intervals based on control neighbours based on bootstrap
#'
#' @param varphi Vector, the integrand on a vector of evaluation points.
#' @param varphi_int Numeric, the integral estimated based on control neighbours.
#' @param eps Numeric, the error critical value.
#' @param s Numeric, the Hölder exponent of the integrand.
#' @param b_out Numeric, number of bootstrap replications to perform.
#' @param d Numeric, the dimension of the design points. Only one dimension is
#' implemented at the moment.
#' @param cdf Function, indicating the cumulative distribution function of the
#' design points.
#' @returns List, containing the following elements:
#' - **$varphi_int** Numeric, the estimated integral value.
#' - **$pi_l** Numeric, the lower bound of the prediction interval.
#' - **$pi_u** Numeric, the upper bound of the prediction interval.
#' - **$width** Numeric, the width of the prediction interval.
#' @export
#'
mc_pi <- function(varphi, varphi_int, eps = 0.05, s,
                  b_out = 200, d = 1, cdf = identity, pow) {

  if(missing(pow)) {
    M_tilde <- floor(length(varphi$t) / 2)
  } else {
    M_tilde <- floor(length(varphi$t)^(pow))
  }


  x_boot <- purrr::map(seq_len(b_out),
                       ~sort(sample.int(n = length(varphi$t),
                                        size = M_tilde,
                                        replace = FALSE)
                             )
                       )

  varphi_boot <- purrr::map(x_boot, ~list(t = varphi$t[.x],
                                          x = varphi$x[.x]))

  varphi_int_boot <- purrr::map_dbl(varphi_boot,
                                    ~mc_int(x = .x$t,
                                            varphi = .x$x,
                                            cdf = cdf)$varphi_int)

  q_boot <- quantile(x = M_tilde^(1/2 + s/d) * (varphi_int_boot - varphi_int),
                     probs = c(eps/2, 1 - eps/2))

  pi_l <- varphi_int + length(varphi$t)^(-1/2 - s/d) * unname(q_boot[1])
  pi_u <- varphi_int + length(varphi$t)^(-1/2 - s/d) * unname(q_boot[2])

  list(varphi_int = varphi_int,
       pi_l = min(pi_l, pi_u),
       pi_u = max(pi_l, pi_u),
       width = abs(pi_u - pi_l))

}


#' Estimate the integral for functional data with bivariate support using
#' control neighbours
#'
#' @param X_list List, containing the following elements:
#' -**$t** Data Frame, with each row containing the coordinates of the sampling
#' point, with columns `t1` and `t2`.
#' -**$x** Vector of observed points, where each point is associated with each
#' row of coordinates in `t`.
#' @returns Vector, with elements containing the estimated integral of each
#' element in `X_list`.
#' @export

mc_int2d <- function(X_list) {

  # Query the 2-NN and obtains ids for LOO-1NN
  nn_id <- purrr::map(X_list, ~RANN::nn2(data = .x$t,
                                         query = .x$t,
                                         k = 2,
                                         eps = 0)$nn.idx[, 2])

  # Compute the control variate term
  d_term <- purrr::map2_dbl(X_list, nn_id, ~mean(.x$x[.y]))

  # Compute Voronoi volumes
  voronoi_list <- purrr::map(X_list, ~deldir::deldir(x = .x$t$t1, y = .x$t$t2))

  tiles_list <- purrr::map(voronoi_list, ~deldir::tile.list(.x))

  vol_list <- purrr::map(tiles_list, function(tile) {
    purrr::map_dbl(tile, ~.x$area)
  })

  # Compute the control variate integral estimate
  v_term <- purrr::map2_dbl(X_list, vol_list, ~sum(.x$x * .y))

  mean_phi <- purrr::map_dbl(X_list, ~mean(.x$x))

  # Compute and return the control neighbour estimate
  mean_phi + v_term - d_term


}








