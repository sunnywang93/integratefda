
#' Constructs an orthonormal cosine basis
#'
#' @param x Vector of evaluation points.
#' @param J Numeric, the number of basis functions to take. Should range
#' between 0 and an integer value, where 0 indicates only the constant function
#' to be taken.
#' @returns Matrix of dimension `x` by `J + 1`.
#' @export
cos_basis <- function(x, J) {
  cbind(rep(1, length(x)),
        sqrt(2) * cos(outer(x, seq_len(J)) * pi))
}

#' Estimate the density using a series expansion
#'
#' @param xi Vector, the observed design points.
#' @param x Vector, the evaluation points to estimate the density.
#' @param J Numeric, the number of basis functions.
#' @returns Vector, the estimated density on the vector of evaluation points.
#' @export

den_series <- function(xi, x, J, thresh = TRUE, cj0 = 3, cj1 = 0.8, cth = 4) {

  if(thresh) {
    Jmax <- floor(cj0 + cj1 * log(length(xi)))
    thetaj <- colMeans(cos_basis(x = xi, J = Jmax))
    var_phi_xi <- apply(cos_basis(x = xi, J = Jmax), 2, var) / length(xi)
    Jmin <- which.min(cumsum(2 * var_phi_xi - thetaj^2))
    phij <- cos_basis(x = x, J = Jmin - 1)
    indic <- thetaj^2 > (cth * var_phi_xi)
    fhat <- c(tcrossprod(thetaj[1:Jmin] * indic[1:Jmin], phij))
  } else {
    phij <- cos_basis(x = x, J = J)
    thetaj <- colMeans(cos_basis(x = xi, J = J))
    fhat <- c(tcrossprod(thetaj, phij))
  }

  fcheck <- fhat - min(min(fhat), -length(xi)^(-1.51))
  if(min(fhat) <= -length(xi)^(-1.51)) {
    fcheck / (thetaj[1] - min(fhat))
  } else {
    fcheck / (thetaj[1] + length(xi)^(-1.51))
  }

}

#' Estimate the mean function of functional data by pooling and smoothing splines
#'
#' @param X_list List, containing the following elements:
#' -**$t** Vector of observed sampling points.
#' -**$x** Vector of observed points.
#' @param x Vector of evaluation points.
#' @returns List, containing the mean function and its evaluation points
#' @export

splines_mean <- function(X_list, x) {

  t_pool <- sort(unlist(purrr::map(X_list, ~.x$t)), index.return = TRUE)
  x_pool <- unlist(purrr::map(X_list, ~.x$x))[t_pool$ix]

  spline_mod <- stats::smooth.spline(x = t_pool$x,
                                     y = x_pool,
                                     cv = FALSE)

  pred_mod <- predict(spline_mod, x = x)

  list(t = pred_mod$x,
       x = pred_mod$y)

}


#' Construct the integrated cosine basis
#'
#' Used for the estimation for the cumulative distribution function.
#'
#' @param x Vector of evaluation points.
#' @param J Numeric, number of basis functions.
#' @returns Matrix, containing the basis functions on the evaluation
#' points.
#' @export

int_cos_basis <- function(x, J) {

  Phi_0 <- 1 - x
  Phi_j <- -sqrt(2) * sin(outer(x, seq_len(J)) * pi) |>
    sweep(MARGIN = 2, STATS = pi * seq_len(J), FUN = "/")

  unname(cbind(Phi_0, Phi_j))

}

#' Estimates the cumulative distribution function using series expansions
#'
#' @param xi Vector of observed points.
#' @param x Vector of evaluation points.
#' @param J Numeric, number of basis functions. Ignored if thresholding
#' is used.
#' @param thresh Boolean, indicating whether thresholding should be performed.
#' @param cj0 Numeric.
#' @param cj1 Numeric.
#' @param cth Numeric.
#' @returns Vector, the estimated cdf.
#' @export

cdf_series <- function(xi, x, J, thresh = TRUE, cj0 = 3, cj1 = 0.8, cth = 4) {

  if(thresh) {
    Jmax <- floor(cj0 + cj1 * log(length(xi)))
    aj <- colMeans(int_cos_basis(x = xi, J = Jmax))
    var_phi_xi <- apply(int_cos_basis(x = xi, J = Jmax), 2, var) / length(xi)
    Jmin <- which.min(cumsum(2 * var_phi_xi - aj^2))
    indic <- aj^2 > (cth * var_phi_xi)
    phij <- cos_basis(x = x, J = Jmin - 1)
    Fhat <- c(tcrossprod(aj[1:Jmin] * indic[1:Jmin], phij))

  } else {
    aj <- colMeans(int_cos_basis(x = xi, J = J))
    Fhat <- c(tcrossprod(aj, cos_basis(x = x, J = J)))
  }

  Fhat

}






