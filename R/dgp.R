#' Generate points from a convex density function
#'
#' Points are generated using the inverse transform sampling method.
#'
#' @param a Numeric, multiplicative parameter.
#' @param n Numeric, number of points to generate
#' @returns Vector of sample points.
#' @export

boat_cdf <- function(a, n) {

  if(a < 0 & a >= 12) {
    stop("a is out of range!")
  }

  b <- 1 - (a / 12)
  U <- runif(n)
  coef_inv <- purrr::map(U, ~c(-.x, a/4 + b, -a/2, a/3))
  roots_inv <- purrr::map(coef_inv, ~polyroot(.x))
  id_real <- purrr::map(roots_inv,
                        ~which((abs(Im(.x)) < 1e-7) & (Re(.x) >= 0) &(Re(.x) < 1)))
  purrr::map2_dbl(roots_inv, id_real, ~Re(.x[.y]))

}
