library(mcscores)
library(parallel)
library(tictoc)
k_length <- 12
t_nvec <- c(50, 100, 200)
nu_vec <- c(1, 1.5, 2)
n_surface <- 500
param_cart <- expand.grid(t_length = t_nvec,
                          nu = nu_vec)
seeds <- sample.int(n = 10000, size = nrow(param_cart))


tic()
result_list <- mclapply(seq_len(nrow(param_cart)), function(prow) {

   t_length <- param_cart[prow, "t_length"]
   nu <- param_cart[prow, "nu"]

   set.seed(seeds[prow])
   # Generate two-dimensional design points
   t_list <- purrr::map(seq_len(n_surface), function(i) {
     t(sapply(seq_len(t_length), function(x) runif(n = 2, 0, 1))) |>
       as.data.frame()
   })

   for(i in 1:length(t_list)) {
     colnames(t_list[[i]]) <- c("t1", "t2")
   }

   # Generate surfaces
   sheets_list <- purrr::map(t_list, ~bs_kl(xout = .x, k = 12, nu1 = nu, nu2 = nu))

   lapply(sheets_list, function(sheet) {
     phi1 <- sqrt(2) * cos(outer(seq_len(k_length), sheet$x_obs$t1) * pi)
     phi2 <- sqrt(2) * cos(outer(seq_len(k_length), sheet$x_obs$t2) * pi)
     phi_big <- sapply(seq_len(ncol(phi1)),
                       function(x) outer(phi1[, x], phi2[, x]),
                       simplify = "array")

     int_sheet <- sheet$x_obs * phi_big[ind_mat[id, 1],ind_mat[id, 2],]
     varphi <- mc_int2d(int_sheet)

   })


   # Construct basis functions
   phi1_list <- purrr::map(sheets_list,
                           ~sqrt(2) * cos(outer(seq_len(k_length),
                                                .x$x_obs$t1) * pi))

   phi2_list <- purrr::map(sheets_list,
                           ~sqrt(2) * cos(outer(seq_len(k_length),
                                                .x$x_obs$t2) * pi))

   ind_mat <- rbind(c(1, 1), c(2, 2), c(3, 3))

   phi_big <- purrr::map2(phi1_list, phi2_list,
                          ~sapply(seq_len(ncol(.x)), function(z) {
                            outer(.x[, z], .y[, z])},
                            simplify = "array"))

   mc_list <- lapply(seq_len(nrow(ind_mat)), function(id) {
     # Construct integrands for scores
     varphi_list <- purrr::map2(sheets_list, phi_big,
                                ~.x$x_obs$x * .y[ind_mat[id, 1],ind_mat[id, 2],]
                                )

     varphi_list2 <- purrr::map2(sheets_list, varphi_list, function(x, y) {
       x$x_obs$x <- y
       return(x)
     })

     purrr::map_dbl(varphi_list2, ~mc_int2d(.x$x_obs))

     # Apply methods and save results
     data.frame(mc_hat = mc_int2d(varphi_list2),
                mu_hat = purrr::map_dbl(varphi_list2, ~mean(.x$x_obs$x)),
                xi_true = purrr::map_dbl(sheets_list,
                                         ~.x$xi[ind_mat[id, 1],ind_mat[id, 2]]),
                group = paste("xi", paste0(ind_mat[id, ], collapse = "."),
                              sep = "_"),
                n = t_length,
                nu = nu
                )
   })

   do.call('rbind', mc_list)

}, mc.cores = detectCores() - 1)
toc()

sum(1 / (outer(seq_len(3)^(2*nu) , seq_len(3)^(2*nu)))) /
  sum(1 / (outer(seq_len(200)^(2*nu) , seq_len(200)^(2*nu))))

# analysis of results ========================================================
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

result_combined <- as_tibble(do.call('rbind', result_list)) |>
  mutate(mc = abs(mc_hat - xi_true) ,
         mu = abs(mu_hat - xi_true))


box_list <- lapply(seq_len(nrow(param_cart)), function(k) {
  result_combined |>
    filter(nu == param_cart[k, "nu"], n == param_cart[k, "t_length"]) |>
    ggplot(aes(x = group, y = ratio)) +
    geom_boxplot() +
    xlab("Scores") +
    ylab("Difference") +
    ggtitle(paste0("M = ",
                   param_cart[k, "t_length"], ", Nu = ",
                   param_cart[k, "nu"])) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
})


mc_ref <- result_combined |>
  filter(n == 200, nu == 1) |>
  select(mc)

mc_1 <- result_combined |>
  filter(n == 200, nu == 1.5) |>
  select(mc)

mc_2 <- result_combined |>
  filter(n == 200, nu == 2) |>
  select(mc)

par(mfrow = c(1, 2))
boxplot(log(mc_1$mc / mc_ref$mc))
boxplot(log(mc_2$mc / mc_ref$mc))

library(patchwork)
wrap_plots(box_list, nrow = 3, ncol = 3) +
  plot_layout(guides = "collect")

wrap_plots(box_list[1:3], nrow = 1, ncol = 3) +
  plot_layout(guides = "collect")
wrap_plots(box_list[4:6], nrow = 1, ncol = 3) +
  plot_layout(guides = "collect")
wrap_plots(box_list[7:9], nrow = 1, ncol = 3) +
  plot_layout(guides = "collect")

