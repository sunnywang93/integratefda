
library(here)
library(dplyr)
library(ggplot2)
library(mcscores)
library(NbClust)
library(tictoc)
# Load data
load(paste0(here(), "/swimmers.rda"))


# Preliminary exploration
head(swimmers)
summary(swimmers)

swimmers |>
  filter(ID %in% paste0(1:10)) |>
  ggplot(aes(x = Age, y = Performance, color = ID)) +
  geom_point() +
  geom_line()

swimmers |>
  group_by(ID) |>
  mutate(n_obs = n()) |>
  summarise(mu_obs = mean(n_obs)) |>
  arrange(desc(mu_obs))

df_norm <- swimmers |>
  mutate(Age = (Age - min(Age)) / (max(Age) - min(Age)))

df_age <- df_norm |>
  select(Age) |>
  arrange(Age)

# estimate the density of the design points
xout <- seq(0, 1, l = 101)

fhat <- list(t = xout,
             x = den_series(xi = df_age$Age,
                            x = xout,
                            thresh = TRUE))

plot(fhat$t, fhat$x, type = "l")

# Map back to normal age
grid_age <- fhat$t * (max(swimmers$Age) - min(swimmers$Age)) + min(swimmers$Age)
plot(grid_age, fhat$x, type = "l")

# Mean estimation
df_list <- df_norm |>
  select(-Gender) |>
  group_by(ID) |>
  group_split(.keep = FALSE) |>
  purrr::map(~list(t = .x$Age,
                   x = .x$Performance))

mu_spl <- splines_mean(X_list = df_list,
                       x = xout)


swim_dense <- df_norm |>
  group_by(ID) |>
  mutate(count = n()) |>
  filter(count >= 25) |>
  select(ID, Age, Performance) |>
  group_split(.keep = FALSE) |>
  purrr::map(~list(t = .x$Age,
                   x = .x$Performance))

tic()
pca_ll <- fdapace::FPCA(Lt = purrr::map(df_list, ~.x$t),
                        Ly = purrr::map(df_list, ~.x$x),
                        optns = list(nRegGrid = 101))
toc()


plot(pca_ll$workGrid, pca_ll$phi[, 1], type = "l")
plot(pca_ll$workGrid, pca_ll$phi[, 2], type = "l")
plot(pca_ll$workGrid, pca_ll$phi[, 3], type = "l")
plot(pca_ll$workGrid, pca_ll$phi[, 4], type = "l")
# PCA with both pooling and basis expansion (with a deterministic basis)
# Inference
# Clustering with the scores, computing using portier vs sample mean

mu_i <- purrr::map(df_list,
                   ~pracma::interp1(x = mu_spl$t,
                                    y = mu_spl$x,
                                    xi = .x$t,
                                    method = "linear"))


phi_i <- purrr::map(df_list,
                    ~apply(pca_ll$phi, 2, function(phij) {
                      pracma::interp1(x = pca_ll$workGrid,
                                      y = phij,
                                      xi = .x$t,
                                      method = "linear")
                    }))

f_i <- purrr::map(df_list, ~pracma::interp1(x = fhat$t,
                                            y = fhat$x,
                                            xi = .x$t,
                                            method = "linear"))


varphi_i <- purrr::map(seq_along(df_list),
                       ~(df_list[[.x]]$x - mu_i[[.x]]) * phi_i[[.x]] / f_i[[.x]]
                       )

cdf_hat <- list(t = xout,
                x = cdf_series(xi = df_age$Age,
                               x = xout,
                               J = 12,
                               thresh = TRUE))




xi_list <- purrr::map2(varphi_i, df_list,
                       ~apply(.x, 2, function(varphi) mc_int(x = .y$t,
                                                             varphi = varphi,
                                                             cdf = cdf_hat)$varphi_int))

df_mc <- do.call('rbind', xi_list)
colnames(df_mc) <- paste0("xi", seq_len(pca_ll$selectK))

xi_emp <- purrr::map2(varphi_i, df_list,
                      ~apply(.x, 2, function(varphi) mean(varphi))
                      )

df_emp <- do.call('rbind', xi_emp)
colnames(df_emp) <- paste0("xi", seq_len(pca_ll$selectK))


scores_df <- as.data.frame(pca_ll$xiEst)
colnames(scores_df) <- paste0("xi", seq_len(pca_ll$selectK))

df_scaled_mc <- scale(df_mc, center = FALSE)
df_scaled_emp <- scale(df_emp, center = FALSE)

df_mc_norm <- sweep(df_mc, 2, sqrt(pca_ll$lambda), FUN = "/")
df_emp_norm <- sweep(df_emp, 2, sqrt(pca_ll$lambda), FUN = "/")

k_opt <- NbClust::NbClust(data = df_mc_norm,
                          method = 'kmeans',
                          index = 'all')

mc_kmeans <- kmeans(x = df_scaled_mc, centers = 2, nstart = 25)
emp_kmeans <- kmeans(x = df_scaled_emp, centers = 2, nstart = 25)
par(mfrow = c(1, 2))
plot(df_scaled_mc, col = mc_kmeans$cluster,
     xlim = c(-10, 15))
plot(df_scaled_emp, col = emp_kmeans$cluster,
     xlim = c(-10, 15))
plot(df_scaled_emp, col = ifelse(emp_kmeans$cluster == 1, 2, 1),
     xlim = c(-10, 15))


mc_kmeans <- kmeans(x = df_mc_norm, centers = 3, nstart = 50)
emp_kmeans <- kmeans(x = df_emp_norm, centers = 3, nstart = 50)
par(mfrow = c(1, 2))
plot(df_mc_norm, col = mc_kmeans$cluster)
plot(df_emp_norm, col = emp_kmeans$cluster)
# try plots after PCA to see in 1D the different clusters

# Extract obs with group 1
df_norm |>
  filter(ID %in% which(mc_kmeans$cluster==1)) |>
  summarise(mu_perf = max(Performance))


df_norm |>
  filter(ID %in% which(mc_kmeans$cluster==2)) |>
  summarise(mu_perf = max(Performance))

?hclust

# CDF can be computed explicitly - used that instead
# Inference - take regularity to be = 1- Fast decrease rate -> Lipschitz

coef_int <- purrr::map(df_list,
                       ~colMeans(int_cos_basis(x = .x$t, J = 12)))

F_hat <- purrr::map(coef_int,
                    ~rowSums(sweep(cos_basis(x = xout, J = 12), 2, .x, FUN = "*")))

