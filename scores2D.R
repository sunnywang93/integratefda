library(mcscores)
library(parallel)
library(tictoc)
library(parallel)
library(doParallel)
library(here)

k_length <- 12
t_nvec <- c(50, 100, 200)
gamma_vec <- c(1, 1.5, 2)
n_surface <- 1000
b_vec <- c(0, 1/2)

param_cart <- expand.grid(t_length = t_nvec,
                          gamma = gamma_vec,
                          b = b_vec)

unif_trans <- function(x, b) (-(1 - b/2) + sqrt((1 - b/2)^2 + 2*b*x)) / b


set.seed(1234)
seeds <- sample.int(n = 10000, size = nrow(param_cart))


n_cores <- detectCores() - 1
cl <- makeCluster(spec = n_cores)
registerDoParallel(cl)

tic()
result_list <- foreach(prow = 1:nrow(param_cart),
                       .packages = 'mcscores') %do% {

   t_length <- param_cart[prow, "t_length"]
   gamma <- param_cart[prow, "gamma"]
   b <- param_cart[prow, "b"]

   set.seed(seeds[prow])
   # Generate two-dimensional design points
   t_list <- purrr::map(seq_len(n_surface), function(i) {
     t(sapply(seq_len(t_length), function(x) runif(n = 2, 0, 1))) |>
       as.data.frame()
   })

   for(i in 1:length(t_list)) {
     colnames(t_list[[i]]) <- c("t1", "t2")
   }

   # Transformation for the non-uniform case - b = 0 corresponds to uniform design
   if(b > 0) {
     t_list <- lapply(seq_along(t_list), function(id) {
       apply(t_list[[id]], 1, function(x) unif_trans(x, b)) |>
         t() |>
         as.data.frame()
     })
   }

   # Generate surfaces
   sheets_list <- purrr::map(t_list,
                             ~bs_kl(xout = .x, k = 12,
                                    gamma1 = gamma, gamma2 = gamma))


   ind_mat <- rbind(c(1, 1), c(2, 2), c(3, 3))

   mc_list <- foreach(s_id = seq_along(sheets_list),
                      .packages = 'mcscores') %dopar% {

     sapply(seq_len(nrow(ind_mat)), function(id) {
       phi1 <- sqrt(2) * cos(outer(seq_len(k_length), sheets_list[[s_id]]$x_obs$t1) * pi)
       phi2 <- sqrt(2) * cos(outer(seq_len(k_length), sheets_list[[s_id]]$x_obs$t2) * pi)
       phi_big <- sapply(seq_len(ncol(phi1)),
                         function(x) outer(phi1[, x], phi2[, x]),
                         simplify = "array")

       # sheets_list[[s_id]]$x_obs$x <- sheets_list[[s_id]]$x_obs$x *
       #   phi_big[ind_mat[id, 1],ind_mat[id, 2],]

       sheets_list[[s_id]]$x_obs$x <- sheets_list[[s_id]]$x_obs$x *
         phi_big[ind_mat[id, 1],ind_mat[id, 2],] /
         (((1 - b/2) + b * sheets_list[[s_id]]$x_obs$t1) *
            ((1 - b/2) + b * sheets_list[[s_id]]$x_obs$t2))

       mc_hat <- mc_int2d(sheets_list[[s_id]]$x_obs)
       mc_pi <- mc_pi2d(X = sheets_list[[s_id]]$x_obs,
                        X_int = mc_hat,
                        eps = 0.05,
                        s = min(gamma - 0.5, 1),
                        b_out = 1000)

       mu_hat <- mean(sheets_list[[s_id]]$x_obs$x, na.rm = TRUE)
       mu_pi <- pi_subsam2d(X = sheets_list[[s_id]]$x_obs,
                            X_int = mu_hat,
                            eps = 0.05,
                            s = min(gamma - 0.5, 1),
                            b_out = 1000,
                            int_fun = mean)

       riemann_hat <- riemann_2d(X = sheets_list[[s_id]]$x_obs)
       riemann_pi <- pi_subsam2d(X = sheets_list[[s_id]]$x_obs,
                                 X_int = riemann_hat,
                                 eps = 0.05,
                                 s = min(gamma - 0.5, 1),
                                 b_out = 1000,
                                 int_fun = riemann_2d)

       c(mc_hat = mc_hat,
         mu_hat = mu_hat,
         riemann_hat = riemann_hat,
         xi_true = sheets_list[[s_id]]$xi_norm[ind_mat[id, 1], ind_mat[id, 2]],
         n = t_length,
         gamma = gamma,
         b = b,
         pi_l = mc_pi$pi_l,
         pi_u = mc_pi$pi_u,
         width = mc_pi$width,
         pi_l_mu = mu_pi$pi_l,
         pi_u_mu = mu_pi$pi_u,
         width_mu = mu_pi$width,
         pi_l_riemann = riemann_pi$pi_l,
         pi_u_riemann = riemann_pi$pi_u,
         width_riemann = riemann_pi$width)
     }) |>
       t()
   }

   xi_names <- paste0("xi",
                      apply(ind_mat, 1, function(x) paste0(x, collapse = ".")))

   result_df <- purrr::map(mc_list, function(x) {
     mc_df <- as.data.frame(x)
     mc_df$group <- xi_names
     mc_df
   })

   result_df

}

toc()
stopCluster(cl)

# analysis of results ========================================================
result_df <- do.call('rbind', purrr::map(result_list, ~do.call('rbind', .x)))

saveRDS(result_df, file = paste0(here(), "/result_df_scores2D.rds"))

sapply(gamma_vec, function(gamma) {
  sum(1 / (outer(seq_len(5)^(2*gamma) , seq_len(5)^(2*gamma)))) /
    sum(1 / (outer(seq_len(200)^(2*gamma) , seq_len(200)^(2*gamma))))
})


library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(xtable)
library(stringr)

result_combined <- as_tibble(result_df) |>
  mutate(mc_err = abs(mc_hat - xi_true),
         mu_err = abs(mu_hat - xi_true),
         riemann_err = abs(riemann_hat - xi_true),
         mu = log(mc_err / mu_err),
         riemann = log(mc_err / riemann_err),
         coverage = (xi_true >= pi_l & xi_true <= pi_u),
         coverage_mu = (xi_true >= pi_l_mu & xi_true <= pi_u_mu),
         coverage_riemann = (xi_true >= pi_l_riemann & xi_true <= pi_u_riemann))


box_list <- lapply(seq_len(nrow(param_cart)), function(k) {
  result_combined |>
  filter(gamma == param_cart[k, "gamma"], n == param_cart[k, "t_length"],
         b == param_cart[k, "b"]) |>
  select(group, mu, riemann) |>
  pivot_longer(cols = c("mu", "riemann"),
               names_to = "Denom",
               values_to = "ratio") |>
  ggplot(aes(x = group, y = ratio, fill = Denom)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, col = "red") +
  scale_fill_grey(start = 0.3, end = 0.9) +
  ylab("Log MAE") +
  xlab("Scores") +
  ggtitle(paste0("M = ",
                 param_cart[k, "t_length"], ", gamma = ",
                 param_cart[k, "gamma"], ", b = ", param_cart[k, "b"])) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
})

# gridExtra::grid.arrange(grobs = box_list)
#
# patchwork::wrap_plots(box_list, ncol = 3, nrow = 3) +
#   patchwork::plot_layout(guides = "collect")
# box_list <- lapply(seq_len(nrow(param_cart)), function(k) {
#   result_combined |>
#     filter(gamma == param_cart[k, "gamma"], n == param_cart[k, "t_length"]) |>
#     ggplot(aes(x = group, y = ratio)) +
#     geom_boxplot() +
#     xlab("Scores") +
#     ylab("Log Absolute Ratio") +
#     ggtitle(paste0("M = ",
#                    param_cart[k, "t_length"], ", gamma = ",
#                    param_cart[k, "gamma"])) +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust = 0.5))
# })


box_titles <- sapply(box_list, function(x) str_replace_all(x$labels[1],
                                                           pattern = " ",
                                                           repl = "") |>
                      str_replace_all(pattern = "=", repl = "") |>
                       str_replace_all(pattern = ",", repl = "_"))

library(tikzDevice)
for(i in 1:length(box_list)) {
  tikz(file = paste0(here(), "/score_box_",
                     box_titles[[i]], ".tex"),
       width = 5,
       height = 5,
       standAlone = TRUE)
  plot(box_list[[i]])
  dev.off()
}


cov_tab <- result_combined |>
  group_by(n, gamma, group, b) |>
  summarise(cov_mc = sum(coverage) / n() * 100,
            cov_mu = sum(coverage_mu) / n() * 100,
            cov_riemann = sum(coverage_riemann) / n() * 100,
            .groups = "keep")

cov_tab$n <- as.integer(cov_tab$n)


width_tab <- result_combined |>
  group_by(n, gamma, group, b) |>
  summarise(width_mc = median(width),
            width_mu = median(width_mu),
            width_riemann = median(width_riemann),
            ratio_mu = (width_mc - width_mu) / width_mu,
            ratio_riemann = (width_mc - width_riemann) / width_riemann,
            .groups = "keep") |>
  select(n, gamma, group, b, ratio_mu, ratio_riemann)


width_tab$n <- as.integer(width_tab$n)


xtable(cov_tab, digits = 1)
xtable(width_tab, digits = 4)
