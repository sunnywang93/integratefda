library(mcscores)
library(parallel)
library(tictoc)
library(reshape2)
library(tidyr)
library(ggplot2)
library(dplyr)

# Parameter settings
k <- 50
xobs_len <- 101
poly_rate <- 3
intercept <- 0
rout <- 2000

# Coefficients for the slope function
phi_coef <- 4 * (-1)^(seq_len(k) + 1) * seq_len(k)^(-poly_rate)

# Confidence intervals without noise ===========================================
conf_level <- 0.95


tic()
result <- mclapply(seq_len(rout), function(rep) {

  # Generate random design points
  xobs <- sort(runif(n = xobs_len))

  # Compute true response using sample path
  xi <- bm_kl_rd(k = k, x = xobs,
                 lambda_rate = poly_rate)
  phi_basis <- sqrt(2) * sin(pi * outer(xobs, seq_len(50) - 0.5))
  yi <- sum(phi_coef * xi$xi)
  slope_obs <- rowSums(sweep(phi_basis, 2, phi_coef, FUN = "*"))

  # Construct integrand
  pdf_list <- list(t = xobs,
                   x = rep(1, length(xobs)))
  varphi <- list(t = xi$t,
                 x = slope_obs * xi$x / pdf_list$x)
  # Estimate response with mean
  yi_hat_mean <- mean(varphi$x)

  # Estimate response with control neighbours
  yi_hat_mc <- mc_int(x = varphi$t,
                      varphi = varphi$x,
                      cdf = identity)

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
  pi_mc <- mc_pi(varphi = varphi,
                 varphi_int = yi_hat_mc$varphi_int,
                 eps = 0.05,
                 s = 0.9,
                 b_out = 200,
                 d = 1,
                 cdf = identity)

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
       yi_hat_mc = yi_hat_mc$varphi_int)

}, mc.cores = 5)
toc()

# Estimation results ===================
mu_diff <- purrr::map_dbl(result, ~.x$mean_diff)
mc_diff <- purrr::map_dbl(result, ~.x$mc_diff)
riemann_diff <- purrr::map_dbl(result, ~.x$riemann_diff)

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



# Confidence intervals with noise ==============================================
conf_level <- 0.95

tic()
result_ci <- mclapply(seq_len(rout), function(rep) {

sigma <- 0.1
xobs <- sort(runif(n = xobs_len))

# Compute true response using sample path
xi <- bm_kl_rd(k = k, x = xobs,
               lambda_rate = poly_rate)

phi_basis <- sqrt(2) * sin(pi * outer(xobs, seq_len(50) - 0.5))
yi <- sum(phi_coef * xi$xi)
slope_obs <- rowSums(sweep(phi_basis, 2, phi_coef, FUN = "*"))

# Construct noisy online set
zi <- list(t = xi$t,
           x = xi$x + sigma * rnorm(n = length(xi$t), mean = 0, sd = 1)
           )
# Estimate response with control neighbours
pdf_list <- list(t = xobs,
                 x = rep(1, length(xobs)))

varphi <- slope_obs * zi$x / pdf_list$x

yi_hat_mc <- mc_int(x = xi$t,
                    varphi = varphi,
                    cdf = identity)

yi_hat_mean <- mean(varphi)

yi_hat_riemann <- pracma::trapz(x = xi$t,
                                y = varphi)

varphi_sigma <- slope_obs * sigma / pdf_list$x

varphi_list <- list(weights = yi_hat_mc$weights,
                    varphi_int = yi_hat_mc$varphi_int,
                    varphi_noise = varphi_sigma)

ci_mc_exact <- mc_ci(varphi_noise = varphi_list,
                     conf_level = conf_level,
                     lim = FALSE)

ci_mc_lim <- mc_ci(varphi_noise = varphi_list,
                   conf_level = conf_level,
                   lim = TRUE)


coverage_exact <- ((yi >= ci_mc_exact$ci_l) & (yi <= ci_mc_exact$ci_u)) |
  ((yi <= ci_mc_exact$ci_l) & (yi >= ci_mc_exact$ci_u))

coverage_lim <- ((yi >= ci_mc_exact$ci_l) & (yi <= ci_mc_exact$ci_u)) |
  ((yi <= ci_mc_exact$ci_l) & (yi >= ci_mc_exact$ci_u))

list(mean_diff = yi_hat_mean - yi,
     mc_diff = yi_hat_mc$varphi_int - yi,
     riemann_diff = yi_hat_riemann - yi,
     ci_exact_count = coverage_exact,
     ci_lim_count = coverage_lim,
     width_ci_exact = ci_mc_exact$width,
     width_ci_lim = ci_mc_lim$width)

}, mc.cores = 5)
toc()

# Estimation results ======
mu_diff <- purrr::map_dbl(result_ci, ~.x$mean_diff)
mc_diff <- purrr::map_dbl(result_ci, ~.x$mc_diff)
riemann_diff <- purrr::map_dbl(result_ci, ~.x$riemann_diff)

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


mu_rel <- purrr::map_dbl(result_ci, ~log(abs(.x$mc_diff / .x$mean_diff)))
rie_rel <- purrr::map_dbl(result_ci, ~log(abs(.x$mc_diff / .x$riemann_diff)))


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



#check coverage
exact_coverage <- sum(purrr::map_lgl(result_ci, ~.x$ci_exact_count))
width_exact <- mean(purrr::map_dbl(result_ci, ~.x$width_ci_exact))

lim_coverage <- sum(purrr::map_lgl(result_ci, ~.x$ci_lim_count))
width_lim <- mean(purrr::map_dbl(result_ci, ~.x$width_ci_lim))



