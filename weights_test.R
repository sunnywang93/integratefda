
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


# x_list <- bm_kl_rd(k = n_basis, x = t_unif)
# y_list <- list(t = t_unif,
#                x = x_list$x,
#                y = x_list$x + rnorm(length(t_unif), 0, 0))


n_basis <- 200
t_unif <- sort(runif(n = 501))
x_list <- bm_kl(k = n_basis, n = 1001)
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

mean(abs(xi_hat$xi_hat - x_list$xi))

# cdf <- list(t = t_unif,
#             x = t_unif)

# cdf_eval <- c(pracma::cumtrapz(x = f_list$t,
#                                y = f_list$x))

n_eval <- 1001
n_random <- 401
n_curves <- 300
points_list <- map(seq_len(n_curves), ~sort(runif(n = n_random)))
n_basis <- 200

x_list <- map(seq_len(n_curves), ~bm_kl(n_basis, n = n_eval))
x_obs_idx <- map(points_list,
                 ~sapply(.x, function(y) which.min(abs(y - seq(0, 1, l = n_eval)))))

y_list <- lapply(seq_along(x_list), function(i) {
  list(t = points_list[[i]],
       x = x_list[[i]]$x[x_obs_idx[[i]]],
       y = x_list[[i]]$x[x_obs_idx[[i]]] + rnorm(length(x_obs_idx[[i]]), 0, 0))
})


mu_list <- map(points_list, ~list(t = .x,
                                  x = rep(0, length(.x))))

psi_list <- map(points_list,
                ~list(t = .x,
                      X = sqrt(2) * sin(outer(pi * .x, seq_len(n_basis) - 0.5)))
)


pdf_list <- map(points_list,
                ~list(t = .x,
                      x = rep(1, length(.x))))

cdf <- identity

xi_hat <- map(seq_len(n_curves), ~mc_scores(y_list = y_list[[.x]],
                                            mu_list = mu_list[[.x]],
                                            psi_list = psi_list[[.x]],
                                            pdf_list = pdf_list[[.x]],
                                            cdf = cdf))

xi_hat_avg <- Reduce('+', map(xi_hat, ~.x$xi_hat)) / length(xi_hat)

xi_true <- Reduce('+', map(x_list, ~.x$xi)) / length(x_list)

mean(abs(xi_hat_avg - xi_true))

