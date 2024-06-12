## Comparison with PACE
sigma <- 0.1
n_basis <- 200
pacecomp <- lapply(seq_len(250), function(rep) {

  t_unif <- sort(runif(n = 101))
  x_list <- bm_kl_rd(k = n_basis, x = t_unif)
  #het_noise <- rnorm(length(t_unif), 0, sd = sigma * (1 + sin(8 * pi * t_unif)/3))
  y_list <- list(t = t_unif,
                 x = x_list$x,
                 y = x_list$x + rnorm(length(t_unif), sd = sigma))

  lambda <- ((seq_len(n_basis) - 0.5) * pi)^(-2)

  mu_list <- list(t = t_unif,
                  x = rep(0, length(t_unif)))


  psi_list <- list(t = t_unif,
                   X = sqrt(2) * sin(outer(pi * t_unif, seq_len(n_basis) - 0.5)))

  Sigma_mat <- lapply(seq_len(n_basis),
                      function(j) lambda[j] * tcrossprod(psi_list$X[, j])) |>
    (\(x) Reduce('+', x))() |>
    (\(x) x + sigma^2 * diag(nrow = nrow(x), ncol = ncol(x)))()

  xi_pace <- pace_scores(y_list, lambda, psi_list, Sigma_mat, mu_list)

  pdf_list <- list(t = t_unif,
                   x = rep(1, length(t_unif)))

  cdf <- identity

  xi_hat <- mc_scores(y_list = y_list,
                      mu_list = mu_list,
                      psi_list = psi_list,
                      pdf_list = pdf_list,
                      cdf = cdf,
                      noise_sd = sigma)

  err_mc <- abs(xi_hat$xi_hat - x_list$xi)
  err_pace <- abs(xi_pace - x_list$xi)

  list(mc_xi = xi_hat$xi_hat,
       pace_xi = xi_pace,
       true_xi = x_list$xi,
       err_mc = err_mc,
       err_pace = err_pace)

})



par(mfrow = c(1, 2))
boxplot(purrr::map_dbl(pacecomp, ~.x$err_mc[1]), main = "MC Score")
boxplot(purrr::map_dbl(pacecomp, ~.x$err_pace[1]), main = "PACE Score")

boxplot(purrr::map_dbl(pacecomp, ~.x$err_mc[2]), main = "MC Score")
boxplot(purrr::map_dbl(pacecomp, ~.x$err_pace[2]), main = "PACE Score")


boxplot(purrr::map_dbl(pacecomp, ~.x$err_mc[3]), main = "MC Score")
boxplot(purrr::map_dbl(pacecomp, ~.x$err_pace[3]), main = "PACE Score")

err_mc <- as.data.frame(t(sapply(pacecomp, function(x) x$err_mc)[1:10, ]))
err_pace <- as.data.frame(t(sapply(pacecomp, function(x) x$err_pace)[1:10, ]))

err_mc_melt <- melt(err_mc)
err_pace_melt <- melt(err_pace)

err_rel <- melt(err_mc / err_pace)

# Plotting the boxplots for the first matrix
plot1 <- ggplot(err_mc_melt, aes(x=variable, y=value)) +
  geom_boxplot() +
  labs(title="Err MC", x="Variable", y="Value")

# Plotting the boxplots for the second matrix
plot2 <- ggplot(err_pace_melt, aes(x=variable, y=value)) +
  geom_boxplot() +
  labs(title="Err PACE", x="Variable", y="Value")

lim_plot <- err_rel |>
  mutate(value_rescaled = quantile(value, c(0, 0.95)))

plot_rel <- ggplot(err_rel, aes(x=variable, y=value)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(err_rel$value, c(0, 0.95))) +
  labs(title="Rel Err", x="Variable", y="Value") +
  geom_hline(yintercept = 1, col = "red")

library(patchwork)

plot1 + plot2


err_mean_mc <- Reduce('+', purrr::map(pacecomp, ~.x$err_mc)) / length(pacecomp)
err_mean_pace <- Reduce('+', purrr::map(pacecomp, ~.x$err_pace)) / length(pacecomp)

