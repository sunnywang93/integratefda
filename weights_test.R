
#' n is like the number of curves to take expectation across
#' m is the number of sampling points
wm_sq_limit <- function(points_list, phi = function(x) 1, cdf = identity) {

  #points_list <- purrr::map(seq_len(n), ~runif(m))

  weights_list <- purrr::map(points_list,
                             ~weights_loo(.x, cdf = cdf))

  weights_fun_prod <- purrr::map2(weights_list, points_list,
                                 ~.x$weights^2 * phi(.y))

  #second_mom <- purrr::map(weights_list, ~.x$weights^2)

  Reduce('+', weights_fun_prod) / length(weights_fun_prod)

}


a <- 1
n <- 100
m <- 25


points_list <- purrr::map(seq_len(1000), function(x) sort(runif(500)))
cdf <- function(x) 1/3*(x - 0.5)^3 + 11*x/12 + 1/24
# some choices of function to play
wm_sq_limit(points_list = points_list,
            phi = function(x) 1,
            cdf = function(x) x)

wm_sq_limit(points_list = points_list,
            phi = identity,
            cdf = cdf) |>
  mean()


wm_sq_limit(points_list = points_list,
            phi = function(x) (sin(x))^2,
            cdf = cdf) |>
  mean()

