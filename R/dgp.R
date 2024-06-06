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
                        ~which((abs(Im(.x)) < 1e-7) & (Re(.x) >= 0) & (Re(.x) < 1)))
  purrr::map2_dbl(roots_inv, id_real, ~Re(.x[.y]))

}


#' Generate brownian motion using the KL decomposition
#'
#' @param k Numeric, number of basis functions.
#' @param n Numeric, number of sampling points.
#' @returns List, containing the sampling and observed points.
#' @export

bm_kl <- function(k, n) {

  x <- seq(0, 1, length.out = n)

  ek <- sqrt(2) * sin(outer(pi * x, seq_len(k) - 0.5))

  lambda_k <- ((seq_len(k) - 0.5) * pi)^(-2)

  xi_k <- rnorm(n = k, mean = 0, sd = 1)

  y <- sweep(ek, 2, sqrt(lambda_k), FUN = "*") |>
    sweep(2, xi_k, FUN = "*") |>
    rowSums()

  list(t = x,
       x = y,
       xi = xi_k)

}


bm_kl_rd <- function(k, x) {

  ek <- sqrt(2) * sin(outer(pi * x, seq_len(k) - 0.5))

  lambda_k <- ((seq_len(k) - 0.5) * pi)^(-2)

  xi_k <- rnorm(n = k, mean = 0, sd = 1)

  y <- sweep(ek, 2, sqrt(lambda_k), FUN = "*") |>
    sweep(2, xi_k, FUN = "*") |>
    rowSums()

  list(t = x,
       x = y,
       xi = xi_k)

}

