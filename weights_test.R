

dir <- "home"

source(paste0(dir, "/weights.R"))


#' n is like the number of curves to take expectation across
#' m is the number of sampling points
wm_sq_limit <- function(n = 100, m = 25, phi = function(x) 1) {

  points_list <- purrr::map(seq_len(n), ~runif(m))

  weights_list <- purrr::map(points_list,
                             ~weights_loo_unif(.x))

  weights_fun_prod <- purrr::map2(weights_list, points_list,
                                 ~.x$weights^2 * phi(.y))

  #second_mom <- purrr::map(weights_list, ~.x$weights^2)

  Reduce('+', weights_fun_prod) / length(weights_fun_prod)

}


# some choices of function to play
wm_sq_limit(n = 1000,
            m = 500,
            phi = function(x) 1)

wm_sq_limit(n = 1000,
            m = 500,
            phi = identity)


wm_sq_limit(n = 1000,
             m = 500,
            phi = function(x) (sin(x))^2)

