
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


points_list <- purrr::map(seq_len(10000), function(x) sort(runif(200)))

purrr::map(points_list, ~degree(.x)$deg^2) |>
  (\(x) Reduce('+', x) / length(x))()

purrr::map(points_list, ~cum_vol(.x)$vol^2) |>
  (\(x) Reduce('+', x) / length(x))()

purrr::map(points_list, ~cum_vol(.x)$vol * degree(.x)$deg)   |>
  (\(x) Reduce('+', x) / length(x))()

purrr::map(points_list, ~cum_vol(.x)$vol^2 + degree(.x)$deg^2 -
             2 * cum_vol(.x)$vol * degree(.x)$deg)   |>
  (\(x) Reduce('+', x) / length(x))()

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




n_basis <- 50
x_list <- bm_kl(k = 50, n = 501)
t_unif <- sort(runif(n = 201))
t_idx <- sapply(t_unif, function(ti) which.min(abs(ti - x_list$t)))
y_list <- list(t = t_unif,
               x = x_list$x[t_idx],
               y = x_list$x[t_idx] + rnorm(length(t_idx), 0, 0))

mu_list <- list(t = t_unif,
                x = rep(0, length(t_unif)))


psi_list <- list(t = t_unif,
                 X = sqrt(2) * sin(outer(pi * t_unif, seq_len(n_basis) - 0.5)))

pdf_list <- list(t = t_unif,
                 x = rep(1, length(t_unif)))

cdf <- identity

xi_hat <- mc_scores(y_list = y_list,
                    mu_list = mu_list,
                    psi_list = psi_list,
                    pdf_list = pdf_list,
                    cdf = cdf)

# cdf <- list(t = t_unif,
#             x = t_unif)

# cdf_eval <- c(pracma::cumtrapz(x = f_list$t,
#                                y = f_list$x))






