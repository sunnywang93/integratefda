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
#' @returns List, containing the sampling points, observed points, and normalised scores.

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
       xi = xi_k * sqrt(lambda_k))

}


#' Generate brownian motion using the KL decomposition on a random grid of points
#'
#' @param k Numeric, number of basis functions.
#' @param x Vector containing the sampling points.
#' @param lambda_rate Numeric, the polynomial rate of decay of eigenvalues.
#' @param norm_constant Numeric, a normalising constant to make the scores centered
#' when using distribution that doesn't have zero mean.
#' @param norm_factor Numeric, a normalising factor to make the scores unit variance
#' when using a distribution that doesn't have unit variance by default.
#' @param xi_dist Function, indicating the distribution in which the scores
#' should be simulated from. Arguments associated to this function should be passed as
#' separate arguments.
#' @returns List, containing the sampling points, observed points, and normalised scores.
#' @export

bm_kl_rd <- function(k, x, lambda_rate, norm_constant, norm_factor, xi_dist, ...) {

  ek <- sqrt(2) * sin(outer(pi * x, seq_len(k) - 0.5))

  lambda_k <- ((seq_len(k) - 0.5) * pi)^(-lambda_rate)

  if(missing(xi_dist)) {
    xi_k <- rnorm(n = k, mean = 0, sd = 1) * sqrt(lambda_k)
  } else {
    xi_k <- (xi_dist(n = k, ...) - norm_constant) * norm_factor * sqrt(lambda_k)
  }

  y <- sweep(ek, 2, xi_k, FUN = "*") |>
    rowSums()

  list(t = x,
       x = y,
       xi = xi_k)

}

#' Simulation brownian sheet using the Kosambi-Karhunen-Loève decomposition
#'
#' @param xout Data frame, consisting of:
#' - **$t1** Vector containing the coordinates in the first dimension.
#' - **$t2** Vector containing the coordinates in the second dimension.
#' @param k Numeric, the number of basis functions.
#' @param gamma1 Numeric, the rate of decay of the square root of eigenvalues in
#' the first dimension.
#' @param gamma2 Numeric, the rate of decay of the square root of eigenvalues in
#' the second dimension.
#' @returns List, containing:
#' - **$x_obs** Data frame containing the observed points at the `(t1, t2)`
#' coordinates.
#' - **$xi** Matrix containing the scores.
#' @export

bs_kl <- function(xout, k, gamma1 = 1, gamma2 = 1) {

  phi1 <- sqrt(2) * cos(outer(seq_len(k), xout$t1) * pi) / (seq_len(k) * pi)^gamma1
  phi2 <- sqrt(2) * cos(outer(seq_len(k), xout$t2) * pi) / (seq_len(k) * pi)^gamma2
  xi <- array(data = rnorm(n = nrow(phi1)^2),
              dim = c(nrow(phi1), nrow(phi1)))

  X_obs <- sapply(seq_len(ncol(phi1)), function(v) {
      crossprod(phi1[, v], xi) %*% phi2[, v]
    })

  list(x_obs = data.frame(t1 = xout$t1,
                          t2 = xout$t2,
                          x = X_obs
                          ),
       xi_norm = xi / outer((seq_len(k) * pi)^gamma1, (seq_len(k) * pi)^gamma2)
       )

}


kl_rd_poly <- function(k, x, lambda_rate) {

  ek <- sqrt(2) * sin(outer(pi * x, seq_len(k) - 0.5))

  lambda_k <- ((seq_len(k) - 0.5) * pi)^(-lambda_rate)

  xi_k1 <- rnorm(n = floor(k/2))
  xi_k2 <- (xi_k1^2 - 1) * sqrt(1/2)
  xi_k <- as.vector(rbind(xi_k1, xi_k2)) * sqrt(lambda_k)

  y <- sweep(ek, 2, xi_k, FUN = "*") |>
    rowSums()

  list(t = x,
       x = y,
       xi = xi_k)

}


