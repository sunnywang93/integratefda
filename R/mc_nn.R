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
      pracma::interp1(x = cdf$t,
                      y = cdf$x,
                      xi = x,
                      method = "linear")
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


#' Computes prediction intervals based on control neighbours based on subsampling
#'
#' @param varphi List, containing the following elements:
#' - **$t** Vector of evaluation points.
#' - **$x** Vector of observed points.
#' @param varphi_int Numeric, the integral estimated based on control neighbours.
#' @param eps Numeric, the error critical value.
#' @param s Numeric, the Hölder exponent of the integrand.
#' @param b_out Numeric, number of bootstrap replications to perform.
#' @param cdf Function or list, indicating the cumulative distribution function of the
#' design points. If list, should contain the following elements:
#' - **$t** Vector of evaluation points.
#' - **$x** Vector of observed points.
#' @returns List, containing the following elements:
#' - **$varphi_int** Numeric, the estimated integral value.
#' - **$pi_l** Numeric, the lower bound of the prediction interval.
#' - **$pi_u** Numeric, the upper bound of the prediction interval.
#' - **$width** Numeric, the width of the prediction interval.
#' @export
#'
mc_pi <- function(varphi, varphi_int, eps = 0.05, s,
                  b_out = 200, cdf = identity) {


  M_tilde <- floor(length(varphi$t) / 2)

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

  q_boot <- quantile(x = M_tilde^(1/2 + s) * (varphi_int_boot - varphi_int),
                     probs = c(eps/2, 1 - eps/2))

  pi_l <- varphi_int + length(varphi$t)^(-1/2 - s) * unname(q_boot[1])
  pi_u <- varphi_int + length(varphi$t)^(-1/2 - s) * unname(q_boot[2])

  list(varphi_int = varphi_int,
       pi_l = min(pi_l, pi_u),
       pi_u = max(pi_l, pi_u),
       width = abs(pi_u - pi_l))

}


#' Estimate the integral for functional data with bivariate support using
#' control neighbours
#'
#' @param X Data frame, containing the following columns:
#' -**$t1** Vector, containing the coordinates in `t1`.
#' -**$t2** Vector, containing the coordinates in `t2`.
#' -**$x** Vector of observed points, where each point is associated with each
#' row of coordinates in `t`.
#' @param xmin Numeric, the minimum coordinate of the `x` axis.
#' @param xmax Numeric, the maximum coordinate of the `x` axis.
#' @param ymin Numeric, the minimum coordinate of the `y` axis.
#' @param ymax Numeric, the maximum coordinate of the `y` axis.
#' @returns Numeric, the estimated integral.
#' @export

mc_int2d <- function(X, xmin = 0, xmax = 1, ymin = 0, ymax = 1) {

  # Query the 2-NN and obtains ids for LOO-1NN
  nn_id <- RANN::nn2(data = X[, c("t1", "t2")],
                     query = X[, c("t1", "t2")],
                     k = 2,
                     eps = 0)$nn.idx[, 2]

  # Compute the control variate term
  d_term <- mean(X[nn_id, ]$x, na.rm = TRUE)

  # Compute Voronoi volumes
  voronoi_vol <- deldir::deldir(x = X$t1,
                                y = X$t2,
                                rw = c(xmin, xmax, ymin, ymax))$summary$dir.area

  # Compute the control variate integral estimate
  v_term <- sum(X$x * voronoi_vol, na.rm = TRUE)

  mean_phi <- mean(X$x, na.rm = TRUE)

  # Compute and return the control neighbour estimate
  mean_phi - d_term + v_term

}


#' Computes prediction intervals for 2D control neighbour integral estimates
#' based on subsampling
#'
#' @param X Data frame, containing the following columns:
#' - **$t1** Vector, containing the coordinates in `t1`.
#' - **$t2** Vector, containing the coordinates in `t2`.
#' - **$x** Vector of observed points, where each point is associated with each
#' row of coordinates in `t`.
#' @param X_int Numeric, the control neighbour integral estimate.
#' @param eps Numeric, the error critical value.
#' @param s Numeric, the Hölder exponent of the integrand.
#' @param b_out Numeric, number of bootstrap replications to perform.
#' @returns List, containing the following elements:
#' - **$X_int** Numeric, the estimated integral value.
#' - **$pi_l** Numeric, the lower bound of the prediction interval.
#' - **$pi_u** Numeric, the upper bound of the prediction interval.
#' - **$width** Numeric, the width of the prediction interval.
#' @export

mc_pi2d <- function(X, X_int, eps = 0.05, s, b_out = 200) {

  M_star <- floor(nrow(X) / 2)

  t_star_id <- sapply(seq_len(b_out), function(b) {
    sample.int(n = nrow(X), size = M_star, replace = FALSE)
  })

  int_boot <- apply(t_star_id, 2, function(b) mc_int2d(X[b, ]))

  q_boot <- quantile(x = M_star^(1/2 + s/2) * (int_boot - X_int),
                     probs = c(eps/2, 1 - eps/2))

  pi_l <- X_int + nrow(X)^(-1/2 - s/2) * unname(q_boot[1])
  pi_u <-  X_int + nrow(X)^(-1/2 - s/2) * unname(q_boot[2])

  list(X_int = X_int,
       pi_l = min(pi_l, pi_u),
       pi_u = max(pi_l, pi_u),
       width = abs(pi_u - pi_l))


}

#' Compute Riemann sum estimates based on Voronoi partition
#'
#' @param X Data frame, containing the following columns:
#' -**$t1** Vector, containing the coordinates in `t1`.
#' -**$t2** Vector, containing the coordinates in `t2`.
#' -**$x** Vector of observed points, where each point is associated with each
#' row of coordinates in `t`.
#' @param xmin Numeric, the minimum coordinate of the `x` axis.
#' @param xmax Numeric, the maximum coordinate of the `x` axis.
#' @param ymin Numeric, the minimum coordinate of the `y` axis.
#' @param ymax Numeric, the maximum coordinate of the `y` axis.
#' @returns Numeric, the estimated integral.
#' @export

riemann_2d <- function(X, xmin = 0, xmax = 1, ymin = 0, ymax = 1) {

  voronoi_vol <- deldir::deldir(x = X$t1,
                                y = X$t2,
                                rw = c(xmin, xmax, ymin, ymax))$summary$dir.area

  sum(X$x * voronoi_vol, na.rm = TRUE)


}

#' Estimates prediction intervals using subsampling
#'
#' Mainly used for comparison purposes. Only works with either the sample mean
#' or riemann sums using the trapezoidal rule.
#'
#' @param varphi Vector, the integrand on a vector of evaluation points.
#' @param varphi_int Numeric, the integral estimated.
#' @param eps Numeric, the error critical value.
#' @param s Numeric, the Hölder exponent of the integrand.
#' @param b_out Numeric, number of bootstrap replications to perform.
#' @param int_fun Function, indicating the method of approximating integrals.
#' Only `mean` and `trapz` function is accepted.
#' @returns List, containing the following elements:
#' - **$varphi_int** Numeric, the estimated integral value.
#' - **$pi_l** Numeric, the lower bound of the prediction interval.
#' - **$pi_u** Numeric, the upper bound of the prediction interval.
#' - **$width** Numeric, the width of the prediction interval.
#' @export
#'

pi_subsam <- function(varphi, varphi_int, eps = 0.05, s,
                      b_out = 200, int_fun) {

  M_tilde <- floor(length(varphi$t) / 2)

  x_boot <- purrr::map(seq_len(b_out),
                       ~sort(sample.int(n = length(varphi$t),
                                        size = M_tilde,
                                        replace = FALSE)
                       )
  )

  varphi_boot <- purrr::map(x_boot, ~list(t = varphi$t[.x],
                                          x = varphi$x[.x]))


  if(identical(int_fun, mean)) {
    varphi_int_boot <- purrr::map_dbl(varphi_boot,
                                      ~mean(x = .x$x))

    q_boot <- quantile(x = M_tilde^(1/2) * (varphi_int_boot - varphi_int),
                       probs = c(eps/2, 1 - eps/2))

    pi_l <- varphi_int + length(varphi$t)^(-1/2) * unname(q_boot[1])
    pi_u <- varphi_int + length(varphi$t)^(-1/2) * unname(q_boot[2])
  } else {
    varphi_int_boot <- purrr::map_dbl(varphi_boot,
                                      ~pracma::trapz(x = .x$t,
                                                     y = .x$x))

    q_boot <- quantile(x = M_tilde^s * (varphi_int_boot - varphi_int),
                       probs = c(eps/2, 1 - eps/2))

    pi_l <- varphi_int + length(varphi$t)^(-s) * unname(q_boot[1])
    pi_u <- varphi_int + length(varphi$t)^(-s) * unname(q_boot[2])
  }


  list(varphi_int = varphi_int,
       pi_l = min(pi_l, pi_u),
       pi_u = max(pi_l, pi_u),
       width = abs(pi_u - pi_l))

}

#' Estimates prediction intervals using subsampling for 2 dimensions.
#'
#' Mainly used for comparison purposes. Only works with either the sample mean
#' or "riemann sums" using voronoi cells.
#'
#' @param X Data frame, containing the following columns:
#' - **$t1** Vector, containing the coordinates in `t1`.
#' - **$t2** Vector, containing the coordinates in `t2`.
#' - **$x** Vector of observed points, where each point is associated with each
#' row of coordinates in `t`.
#' @param X_int Numeric, the control neighbour integral estimate.
#' @param eps Numeric, the error critical value.
#' @param s Numeric, the Hölder exponent of the integrand.
#' @param b_out Numeric, number of bootstrap replications to perform.
#' @param int_fun Function, indicating the method of approximating integrals.
#' Only `mean` and `riemann_2d` function is accepted.
#' @returns List, containing the following elements:
#' - **$X_int** Numeric, the estimated integral value.
#' - **$pi_l** Numeric, the lower bound of the prediction interval.
#' - **$pi_u** Numeric, the upper bound of the prediction interval.
#' - **$width** Numeric, the width of the prediction interval.
#' @export

pi_subsam2d <- function(X, X_int, eps = 0.05, s, b_out = 200, int_fun) {

  M_star <- floor(nrow(X) / 2)

  t_star_id <- sapply(seq_len(b_out), function(b) {
    sample.int(n = nrow(X), size = M_star, replace = FALSE)
  })


  if(identical(int_fun, mean)) {
    int_boot <- apply(t_star_id, 2, function(b) mean(X[b, ]$x))
    q_boot <- quantile(x = M_star^(1/2) * (int_boot - X_int),
                       probs = c(eps/2, 1 - eps/2))

    pi_l <- X_int + nrow(X)^(-1/2) * unname(q_boot[1])
    pi_u <-  X_int + nrow(X)^(-1/2) * unname(q_boot[2])
  } else {
    int_boot <- apply(t_star_id, 2, function(b) riemann_2d(X[b, ]))
    q_boot <- quantile(x = M_star^s * (int_boot - X_int),
                       probs = c(eps/2, 1 - eps/2))

    pi_l <- X_int + nrow(X)^(-s) * unname(q_boot[1])
    pi_u <-  X_int + nrow(X)^(-s) * unname(q_boot[2])
  }

  list(X_int = X_int,
       pi_l = min(pi_l, pi_u),
       pi_u = max(pi_l, pi_u),
       width = abs(pi_u - pi_l))
}





