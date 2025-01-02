library(integratefda)
library(parallel)
library(doParallel)
library(tictoc)
library(reshape2)
library(tidyr)
library(ggplot2)
library(dplyr)
library(stringr)
library(xtable)

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
conf_level <- 0.95


set.seed(1234)
seeds <- sample.int(n = 10000, size = rout)


n_cores <- detectCores() - 1
cl <- makeCluster(spec = n_cores)
registerDoParallel(cl)


tic()
result_ci <- foreach(prow = 1:nrow(param_cart),
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

  sigma <- 0.1
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

  # Construct noisy online set
  zi <- list(t = xi$t,
             x = xi$x + sigma * rnorm(n = length(xi$t), mean = 0, sd = 1)
  )
  # Construct integrand
  pdf_list <- list(t = xobs,
                   x = 1 - b/2 + b*xobs)

  varphi <- list(t = zi$t,
                 x = slope_obs * zi$x / pdf_list$x)

  # Estimation
  if(b > 0) {
    yi_hat_mc <- mc_int(x = varphi$t,
                        varphi = varphi$x,
                        cdf = function(x) (1 - b/2)*x + b*x^2/2)
  } else {
    yi_hat_mc <- mc_int(x = varphi$t,
                        varphi = varphi$x,
                        cdf = identity)
  }

  yi_hat_mean <- mean(varphi$x)

  yi_hat_riemann <- pracma::trapz(x = varphi$t,
                                  y = varphi$x * pdf_list$x)

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
       width_ci_lim = ci_mc_lim$width,
       n = xobs_len,
       nu = poly_rate,
       b = b)
  }
}
toc()
stopCluster(cl)

# Estimation results ======
result_df <- do.call('rbind', purrr::map(result_ci, ~do.call('rbind', .x))) |>
  apply(2, unlist) |>
  as_tibble() |>
  mutate(mu_rel = log(abs(mc_diff / mean_diff)),
         rie_rel = log(abs(mc_diff / riemann_diff)))


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
    rename(Mean = mu_rel, Riemann = rie_rel) |>
    pivot_longer(cols = Mean:Riemann,
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
  tikz(file = paste0(here(), "/reg_box_noisy_",
                     box_titles[[i]], ".tex"),
       width = 5,
       height = 5,
       standAlone = TRUE)
  plot(box_logmae_list[[i]])
  dev.off()
}


cov_tab <- result_df |>
  select(n, nu, b,
         ci_exact_count, ci_lim_count) |>
  rename(exact = ci_exact_count,
         lim = ci_lim_count) |>
  group_by(n, nu, b) |>
  summarise(coverage_exact = sum(exact) / rout,
            coverage_lim = sum(lim) / rout)


width_tab <- result_df |>
  select(n, nu, b,
         width_ci_exact, width_ci_lim) |>
  rename(exact = width_ci_exact,
         lim = width_ci_lim) |>
  group_by(n, nu, b) |>
  summarise(width_exact = mean(exact),
            width_lim = mean(lim))

left_join(cov_tab, width_tab, by = c("n", "nu", "b")) |>
  mutate(n = as.integer(n),
         nu = as.integer(nu)) |>
  (\(x) print(xtable(x, digits = 3), include.rownames=FALSE))()



