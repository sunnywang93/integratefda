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



