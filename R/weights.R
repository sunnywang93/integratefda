#' Calculates the degree of a vector
#'
#' For each point in a vector, the number of points in which it is
#' a nearest neighbour to is calculated.
#'
#' @param x Vector of points.
#' @returns List, containing the associated degree for each sampling
#' point.
#' @export


degree <- function(x) {

  x_sort <- sort(x, index.return = TRUE)

  dx1 <- diff(x_sort$x[-1])
  dx2 <- diff(x_sort$x)[-(length(x_sort$x) - 1)]

  boo1 <- (dx1 <= dx2) * 1
  boo2 <- 1 - boo1

  dm1 <- boo1[1:(length(boo1)-2)]
  dm2 <- boo2[3:length(boo2)]

  d1 <- (dx2[1] <= dx1[1]) * 1
  dlast <- (dx1[length(dx1)] <= dx2[length(dx2)]) * 1

  d2 <- ((dx1[1] <= dx1[2]) * 1) + 1
  dpen <- ((dx1[length(dx1) - 1] <= dx1[length(dx1) - 2]) * 1) + 1

  deg <- c(d1, d2, dm1 + dm2, dpen, dlast)

  list(x = x,
       deg = deg[order(x_sort$ix)])


}

#' Calculates the leave-one-out cumulative volume for general distributions
#'
#' Given the distribution function and a vector of points, the leave-one-out
#' cumulative volume in one dimension is calculated based on order statistics.
#'
#' @param x Vector of points.
#' @param cdf Function, representation the empirical distribution function.
#' @returns List, containing the associated leave-one-out cumulative volume
#' for each sampling point.
#' @export

cum_vol <- function(x, cdf = identity) {

  x_sort <- sort(x, index.return = TRUE)

  rollsum <- RcppRoll::roll_sum(x = x_sort$x, n = 2)

  id_int <- seq(3, length(x) - 2)

  cm1 <- cdf(rollsum[id_int] / 2) - cdf(rollsum[seq(2, length(x) - 3)] / 2)

  id_plus2 <- id_int + 2
  id_minus2 <- id_int - 2

  cm2 <- cdf(rollsum[id_int] / 2) -
    cdf((x_sort$x[id_int] + x_sort$x[id_minus2]) / 2)
  cm3 <- cdf((x_sort$x[id_int] + x_sort$x[id_plus2]) / 2) -
    cdf(rollsum[seq(2, length(x) - 3)] / 2)

  cm <- (length(x) - 3) * cm1 + cm2 + cm3

  c1 <- (length(x) - 2) * cdf((x_sort$x[1] + x_sort$x[2]) / 2) +
    cdf((x_sort$x[1] + x_sort$x[3]) / 2)

  c21 <- cdf((x_sort$x[2] + x_sort$x[3]) / 2) -
    cdf((x_sort$x[1] + x_sort$x[2]) / 2)
  c2 <- (length(x) - 2) * c21 + cdf((x_sort$x[2] + x_sort$x[4]) / 2)

  cpen1 <- cdf((x_sort$x[length(x) - 1] + x_sort$x[length(x)]) / 2) -
    cdf((x_sort$x[length(x) - 2] + x_sort$x[length(x) - 1]) / 2)
  cpen <- (length(x) - 2) * cpen1 +
    (1 - cdf((x_sort$x[length(x) - 3] + x_sort$x[length(x - 1)]) / 2))

  clast1 <- 1 - cdf((x_sort$x[length(x) - 1] + x_sort$x[length(x)]) / 2)
  clast2 <- 1 - cdf((x_sort$x[length(x) - 2] + x_sort$x[length(x)]) / 2)
  clast <- (length(x) - 2) * clast1 + clast2

  vol <- c(c1, c2, cm, cpen, clast)

  list(x = x,
       vol = vol[order(x_sort$ix)])


}


#' Calculates the leave-one-out cumulative volume for the uniform distribution
#'
#' @param x Vector of points.
#' @returns List, containing the points and its associated cumulative
#' volume.
#' @export

cum_vol_unif <- function(x) {

  x_sort <- sort(x, index.return = TRUE)

  cm1 <- diff(x_sort$x, lag = 2) / 2
  cm2 <- diff(x_sort$x, lag = 4) / 2
  cm <- (length(x) - 2) * cm1[2:(length(cm1)-1)] +
    cm2

  c1 <- (x_sort$x[2] + x_sort$x[1]) * (length(x) - 2) / 2 +
    (x_sort$x[1] + x_sort$x[3]) / 2
  clast1 <- 1 - (x_sort$x[length(x) - 1] + x_sort$x[length(x)]) / 2
  clast2 <- 1 - (x_sort$x[length(x) - 2] + x_sort$x[length(x)]) / 2
  clast <- (length(x) - 2) * clast1 + clast2

  c2 <- (x_sort$x[3] - x_sort$x[1]) * (length(x) - 2) / 2 +
    (x_sort$x[2] + x_sort$x[4]) / 2

  cpen1 <- (x_sort$x[length(x)] - x_sort$x[length(x) - 2]) *
    (length(x) - 2) / 2
  cpen2 <- 1 - (x_sort$x[length(x) - 1] + x_sort$x[length(x) - 3]) / 2
  cpen <-  cpen1 + cpen2

  vol <- c(c1, c2, cm, cpen, clast)

  list(x = x,
       vol = vol[order(x_sort$ix)])

}

#' Compute weights for leave-one-out control neighbours
#'
#' @param x Vector of sampling points.
#' @param cdf Function, representing the cumulative distribution function.
#' @returns List, containing the associated weights for each sampling point.
#' @export
#'
weights_loo <- function(x, cdf = identity) {

  cm_list <- cum_vol(x, cdf = cdf)
  dm_list <- degree(x)

  wm <- (1 + cm_list$vol - dm_list$deg) / length(x)

  list(x = x,
       weights = wm)

}


#' Calculates the leave-one-out-weights for the uniform distribution
#'
#' @param x Vector of points.
#' @returns List, containing the points and asssociated weights.
#' @export

weights_loo_unif <- function(x) {

  cm_list <- cum_vol_unif(x)
  dm_list <- degree(x)

  wm <- (1 + cm_list$vol - dm_list$deg) / length(x)

  list(x = x,
       weights = wm)

}




