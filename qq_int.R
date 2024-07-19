xobs_len <- 101
k <- 50
poly_rate <- 3
intercept <- 0
rout <- 2000
phi_coef <- 4 * (-1)^(seq_len(k) + 1) * seq_len(k)^(-poly_rate)

# Compute true response using sample path
x_dense = seq(0, 1, l = 10000)
xi <- bm_kl_rd(k = k, x = x_dense,
               lambda_rate = poly_rate)
phi_basis <- sqrt(2) * sin(pi * outer(x_dense, seq_len(50) - 0.5))
yi <- sum(phi_coef * xi$xi)
slope_dense <- rowSums(sweep(phi_basis, 2, phi_coef, FUN = "*"))



test <- lapply(seq_len(rout), function(i) {


xobs <- sort(runif(n = xobs_len))

# Construct integrand
pdf_list <- list(t = xobs,
                 x = rep(1, length(xobs)))

slope_obs <- pracma::interp1(x = x_dense,
                             y = slope_dense,
                             xi = xobs)

xi_obs <- pracma::interp1(x = x_dense,
                          y = xi$x,
                          xi = xobs)

varphi <- list(t = xobs,
               x = slope_obs * xi_obs / pdf_list$x)

# Estimate response with control neighbours
yi_hat_mc <- mc_int(x = xobs,
                    varphi = varphi$x,
                    cdf = identity)


list(int_test = yi_hat_mc$varphi_int,
     int_cen = yi_hat_mc$varphi_int - yi,
     int_stand = length(xobs)^(1/2 + (1 - log(length(xobs))^(-2)))*
       (yi_hat_mc$varphi_int - yi))

})



test_vec <- purrr::map_dbl(test, ~.x$int_test)
test_vec2 <- purrr::map_dbl(test, ~.x$int_cen)
test_vec3 <- purrr::map_dbl(test, ~.x$int_stand)

qqnorm(test_vec)
qqline(test_vec)

qqnorm(test_vec2)
qqline(test_vec2)

qqnorm(test_vec3)
qqline(test_vec3)

