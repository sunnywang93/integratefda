library(mcscores)
library(parallel)
library(doParallel)
library(tictoc)
library(reshape2)
library(tidyr)
library(ggplot2)
library(dplyr)

# Parameter settings
k <- 50
xobs_vec <- c(50, 100, 200)
poly_vec <- c(2, 3, 4)
b_vec <- c(0, 0.5)

param_cart <- expand.grid(xobs_len = xobs_vec,
                          poly_rate = poly_vec,
                          b = b_vec)
intercept <- 0
rout <- 500
unif_trans <- function(x, b) (-(1 - b/2) + sqrt((1 - b/2)^2 + 2*b*x)) / b

# Confidence intervals without noise ===========================================
conf_level <- 0.95

set.seed(1234)
seeds <- sample.int(n = 10000, size = rout)


n_cores <- detectCores() - 1
cl <- makeCluster(spec = n_cores)
registerDoParallel(cl)

tic()
result <- foreach(prow = 1:nrow(param_cart),
                  .packages = "mcscores") %do% {


  # set varying parameters
  xobs_len <- param_cart[prow, "xobs_len"]
  poly_rate <- param_cart[prow, "poly_rate"]
  b <- param_cart[prow, "b"]

  # Coefficients for the slope function
  phi_coef <- 4 * (-1)^(seq_len(k) + 1) * seq_len(k)^(-poly_rate)

  # Generate random design points
  result_list <- foreach(rep = 1:rout, .packages = "mcscores") %dopar% {

    set.seed(seeds[rep])

    # Generate design points
    xobs <- sort(runif(n = xobs_len))

    if(b > 0) {
      xobs <- unif_trans(x = xobs, b = b)
    }

    # Compute true response using sample path
    xi <- bm_kl_rd(k = k, x = xobs,
                   lambda_rate = poly_rate)
    phi_basis <- sqrt(2) * sin(pi * outer(xobs, seq_len(50) - 0.5))
    yi <- sum(phi_coef * xi$xi)
    slope_obs <- rowSums(sweep(phi_basis, 2, phi_coef, FUN = "*"))

    # Construct integrand
    pdf_list <- list(t = xobs,
                     x = 1 - b/2 + b*xobs)

    varphi <- list(t = xi$t,
                   x = slope_obs * xi$x / pdf_list$x)
    # Estimate response with mean
    yi_hat_mean <- mean(varphi$x)

    # Estimate response with control neighbours
    if(b > 0) {
      yi_hat_mc <- mc_int(x = varphi$t,
                          varphi = varphi$x,
                          cdf = function(x) (1 - b/2)*x + b*x^2/2)
    } else {
      yi_hat_mc <- mc_int(x = varphi$t,
                          varphi = varphi$x,
                          cdf = identity)
    }

    # Estimate response with Riemann sums
    yi_hat_riemann <- pracma::trapz(x = varphi$t,
                                    y = varphi$x)

      # Construct prediction intervals for mean
    sm <- length(varphi$x)^(-1/2) * sqrt(var(varphi$x))
    z_value <- qnorm(1 - (1 - conf_level)/2)
    ci_l_mean <- yi_hat_mean - z_value * sm
    ci_u_mean <- yi_hat_mean + z_value * sm
    width_ci_mean <- abs(ci_u_mean - ci_l_mean)
    coverage_mean <- ((yi >= ci_l_mean) & (yi <= ci_u_mean)) |
      ((yi <= ci_l_mean) & (yi >= ci_u_mean))

    # Construct prediction intervals for control neighbours
    if(b > 0) {
      pi_mc <- mc_pi(varphi = varphi,
                     varphi_int = yi_hat_mc$varphi_int,
                     eps = 0.05,
                     s = min((poly_rate - 1) / 2, 1),
                     b_out = 200,
                     cdf = function(x) (1 - b/2)*x + b*x^2/2)
    } else {
      pi_mc <- mc_pi(varphi = varphi,
                     varphi_int = yi_hat_mc$varphi_int,
                     eps = 0.05,
                     s = min((poly_rate - 1) / 2, 1),
                     b_out = 200,
                     cdf = identity)
    }

    coverage_mc <- ((yi >= pi_mc$pi_l) & (yi <= pi_mc$pi_u)) |
      ((yi <= pi_mc$pi_l) & (yi >= pi_mc$pi_u))


    list(mean_diff = yi_hat_mean - yi,
         mc_diff = yi_hat_mc$varphi_int - yi,
         riemann_diff = yi_hat_riemann - yi,
         coverage_mean = coverage_mean,
         coverage_mc = coverage_mc,
         yi_hat_mean = yi_hat_mean,
         width_ci_mean = width_ci_mean,
         width_mc = pi_mc$width,
         yi_hat_mc = yi_hat_mc$varphi_int,
         n = xobs_len,
         nu = poly_rate,
         b = b)
  }
  result_list
}
toc()
stopCluster(cl)

# Estimation results ===================
mu_diff <- purrr::map_dbl(result[[18]], ~.x$mean_diff)
mc_diff <- purrr::map_dbl(result[[18]], ~.x$mc_diff)
riemann_diff <- purrr::map_dbl(result[[18]], ~.x$riemann_diff)

result_df <- do.call('rbind', purrr::map(result, ~do.call('rbind', .x))) |>
  apply(2, unlist) |>
  as_tibble() |>



df_diff <- data.frame(
  values = c(mu_diff, mc_diff, riemann_diff),
  group = factor(rep(c("Mean", "NN", "Riemann"), each = rout))
)

# Sanity check
all.equal(filter(df_diff, group == "Mean")$values, mu_diff)
all.equal(filter(df_diff, group == "NN")$values, mc_diff)
all.equal(filter(df_diff, group == "Riemann")$values, riemann_diff)

ggplot(df_diff, aes(x = group, y = values, fill = group)) +
  geom_boxplot() +
  labs(title = "Bias-Variance Comparison", x = "Method", y = "Difference") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_grey(start = 0.3, end = 0.99)



mu_rel <- purrr::map_dbl(result, ~log(abs(.x$mc_diff / .x$mean_diff)))
rie_rel <- purrr::map_dbl(result, ~log(abs(.x$mc_diff / .x$riemann_diff)))


df <- data.frame(
  values = c(mu_rel, rie_rel),
  Group = factor(rep(c("Mean", "Riemann"), each = rout))
)


# Create the boxplots using ggplot2 with different colors
ggplot(df, aes(x = Group, y = values, fill = Group)) +
  geom_boxplot() +
  labs(title = "Mean Absolute Ratio (Log scale)", x = "Method", y = "Err") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_grey(start = 0.4, end = 0.9)


# Inference results ===================
mc_cov <- purrr::map_lgl(result, ~.x$coverage_mc)
width_mc <- purrr::map_dbl(result, ~.x$width_mc)
width_mean <- purrr::map_dbl(result, ~.x$width_ci_mean)
mean_cov <- purrr::map_dbl(result, ~.x$coverage_mean)

#Check coverage and average widths
sum(mc_cov) / rout
sum(mean_cov) / rout
mean(width_mc)
mean(width_mean)


