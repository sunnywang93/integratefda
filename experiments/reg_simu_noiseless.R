library(integratefda)
library(parallel)
library(doParallel)
library(tictoc)
library(reshape2)
library(tidyr)
library(ggplot2)
library(dplyr)
library(stringr)

# Parameter settings
k <- 50
xobs_vec <- c(50, 100, 200)
poly_vec <- c(2, 3, 4)
b_vec <- c(0, 0.5)

param_cart <- expand.grid(xobs_len = xobs_vec,
                          poly_rate = poly_vec,
                          b = b_vec)
intercept <- 0
rout <- 2000
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
                  .packages = "integratefda") %do% {


  # set varying parameters
  xobs_len <- param_cart[prow, "xobs_len"]
  poly_rate <- param_cart[prow, "poly_rate"]
  b <- param_cart[prow, "b"]

  # Coefficients for the slope function
  phi_coef <- 4 * (-1)^(seq_len(k) + 1) * seq_len(k)^(-poly_rate)

  # Generate random design points
  result_list <- foreach(rep = 1:rout, .packages = "integratefda") %dopar% {

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

    # Estimate response with Riemann sums - we need to multiply back
    # the integrnad because Riemann sums does not require to divide by the
    # density since it already converges to the integral
    yi_hat_riemann <- pracma::trapz(x = varphi$t,
                                    y = varphi$x * pdf_list$x)

      # Construct prediction intervals for mean
    sm <- length(varphi$x)^(-1/2) * sqrt(var(varphi$x))
    z_value <- qnorm(1 - (1 - conf_level)/2)
    ci_l_mean <- yi_hat_mean - z_value * sm
    ci_u_mean <- yi_hat_mean + z_value * sm
    width_ci_mean <- abs(ci_u_mean - ci_l_mean)
    coverage_mean <- ((yi >= ci_l_mean) & (yi <= ci_u_mean))

    # Construct prediction intervals for control neighbours
    if(b > 0) {
      pi_mc <- mc_pi(varphi = varphi,
                     varphi_int = yi_hat_mc$varphi_int,
                     eps = 0.05,
                     s = min((poly_rate - 1) / 2, 1),
                     b_out = 1000,
                     cdf = function(x) (1 - b/2)*x + b*x^2/2)
    } else {
      pi_mc <- mc_pi(varphi = varphi,
                     varphi_int = yi_hat_mc$varphi_int,
                     eps = 0.05,
                     s = min((poly_rate - 1) / 2, 1),
                     b_out = 1000,
                     cdf = identity)
    }

    # Prediction intervals for trapezoidal rule based on subsampling
    pi_riemann <- pi_subsam(varphi = list(t = varphi$t,
                                          x = varphi$x * pdf_list$x),
                            varphi_int = yi_hat_riemann,
                            eps = 0.05,
                            s = min((poly_rate - 1) / 2, 1),
                            b_out = 1000,
                            int_fun = pracma::trapz)

    # Prediction intervals for mean based on subsampling
    pi_mean <- pi_subsam(varphi = varphi,
                         varphi_int = yi_hat_mean,
                         eps = 0.05,
                         s = min((poly_rate - 1) / 2, 1),
                         b_out = 1000,
                         int_fun = mean)

    coverage_mc <- ((yi >= pi_mc$pi_l) & (yi <= pi_mc$pi_u))

    coverage_riemann <- ((yi >= pi_riemann$pi_l) & (yi <= pi_riemann$pi_u))

    coverage_mu_sub <- ((yi >= pi_mean$pi_l) & (yi <= pi_mean$pi_u))


    list(mean_diff = yi_hat_mean - yi,
         mc_diff = yi_hat_mc$varphi_int - yi,
         riemann_diff = yi_hat_riemann - yi,
         coverage_mean = coverage_mean,
         coverage_mc = coverage_mc,
         coverage_riemann = coverage_riemann,
         coverage_mu_sub = coverage_mu_sub,
         width_ci_mean = width_ci_mean,
         width_mc = pi_mc$width,
         width_riemann = pi_riemann$width,
         width_mc_sub = pi_mean$width,
         n = xobs_len,
         nu = poly_rate,
         b = b)
  }
  result_list
}
toc()
stopCluster(cl)

# Estimation results ===================
result_df <- do.call('rbind', purrr::map(result, ~do.call('rbind', .x))) |>
  apply(2, unlist) |>
  as_tibble() |>
  mutate(mu_rel = log(abs(mc_diff / mean_diff)),
         rie_rel = log(abs(mc_diff / riemann_diff)))

saveRDS(result_df, file = paste0(here(), "/result_df_reg1D_noiseless.rds"))
result_df <- readRDS(file = paste0(here(), "/result_df_reg1D_noiseless.rds"))

box_diff_list <- lapply(seq_len(nrow(param_cart)), function(k) {
  result_df |>
    filter(n == param_cart[k, "xobs_len"],
           nu == param_cart[k, "poly_rate"],
           b == param_cart[k, "b"]) |>
    select(mean_diff, mc_diff, riemann_diff) |>
    rename(mean = mean_diff, mc = mc_diff, riemann = riemann_diff) |>
    pivot_longer(cols = mean:riemann,
                 names_to = "Method",
                 values_to = "Diff") |>
    ggplot(aes(x = Method, y = Diff, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, col = "red") +
    ggtitle(paste0("M = ",
                   param_cart[k, "xobs_len"], ", nu = ",
                   param_cart[k, "poly_rate"], ", b = ", param_cart[k, "b"])) +
    xlab("Method") +
    ylab("Difference") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_grey(start = 0.3, end = 0.99)
})


box_logmae_list <- lapply(seq_len(nrow(param_cart)), function(k) {
  result_df |>
    filter(n == param_cart[k, "xobs_len"],
           nu == param_cart[k, "poly_rate"],
           b == param_cart[k, "b"]) |>
    select(mu_rel, rie_rel) |>
    rename(mean = mu_rel, riemann = rie_rel) |>
    pivot_longer(cols = mean:riemann,
                 names_to = "Method",
                 values_to = "Log_MAE") |>
    ggplot(aes(x = Method, y = Log_MAE, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, col = "red") +
    ggtitle(paste0("M = ",
                   param_cart[k, "xobs_len"], ", nu = ",
                   param_cart[k, "poly_rate"], ", b = ", param_cart[k, "b"])) +
    xlab("Method") +
    ylab("Log MAE") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_grey(start = 0.3, end = 0.99)
})

box_titles <- sapply(box_logmae_list, function(x) str_replace_all(x$labels$title,
                                                                  pattern = " ",
                                                                  repl = "") |>
                       str_replace_all(pattern = "=", repl = "") |>
                       str_replace_all(pattern = ",", repl = "_"))

library(tikzDevice)
for(i in 1:length(box_logmae_list)) {
  tikz(file = paste0(here(), "/reg_box_",
                     box_titles[[i]], ".tex"),
       width = 5,
       height = 5,
       standAlone = TRUE)
  plot(box_logmae_list[[i]])
  dev.off()
}

# df <- data.frame(
#   values = c(mu_rel, rie_rel),
#   Group = factor(rep(c("Mean", "Riemann"), each = rout))
# )
#
#
# # Create the boxplots using ggplot2 with different colors
# ggplot(df, aes(x = Group, y = values, fill = Group)) +
#   geom_boxplot() +
#   labs(title = "Mean Absolute Ratio (Log scale)", x = "Method", y = "Err") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_grey(start = 0.4, end = 0.9)


# Inference results ===================
cov_tab <- result_df |>
  select(n, nu, b,
         coverage_mean, coverage_mc, coverage_riemann, coverage_mu_sub) |>
  rename(mean = coverage_mean,
         mc = coverage_mc,
         riemann = coverage_riemann,
         mean_sub = coverage_mu_sub) |>
  group_by(n, nu, b) |>
  summarise(coverage_mean = sum(mean) / rout,
            coverage_mc = sum(mc) / rout,
            coverage_riemann = sum(riemann) / rout,
            coverage_mu_sub = sum(mean_sub) / rout)


width_tab <- result_df |>
  select(n, nu, b,
         width_ci_mean, width_mc, width_riemann, width_mc_sub) |>
  rename(mu = width_ci_mean,
         mc = width_mc,
         riemann = width_riemann,
         mean_sub = width_mc_sub) |>
  group_by(n, nu, b) |>
  summarise(width_mean = mean(mu),
            width_mc = mean(mc),
            width_riemann = mean(riemann),
            width_mu_sub = mean(mean_sub))


left_join(cov_tab, width_tab, by = c("n", "nu", "b")) |>
  mutate(n = as.integer(n),
         nu = as.integer(nu)) |>
  (\(x) print(xtable(x, digits = 2), include.rownames=FALSE))()


