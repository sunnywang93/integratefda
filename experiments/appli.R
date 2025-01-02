
library(here)
library(dplyr)
library(ggplot2)
library(integratefda)
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
  geom_line() +
  ggtitle('Performance curves for 10 swimmers') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

swimmers |>
  group_by(ID) |>
  mutate(n_obs = n()) |>
  select(n_obs) |>
  arrange(desc(n_obs)) |>
  distinct()


swimmers |>
  group_by(ID) |>
  mutate(n_obs = n()) |>
  ungroup() |>
  summarise(mu_obs = mean(n_obs))

df_norm <- swimmers |>
  mutate(Age = (Age - min(Age)) / (max(Age) - min(Age))) |>
  mutate(ID = as.numeric(ID)) |>
  as_tibble()

df_age <- df_norm |>
  select(Age) |>
  arrange(Age)


# =============================================================================

# set the grid to be 1/120 - corresponding to identifying age by the month over
# 10 years (120 months)
xout <- seq(0, 1, by = 1/120)

# density estimation by efromovich
fhat <- list(t = xout,
             x = den_series(xi = df_age$Age,
                            x = xout,
                            thresh = TRUE))

# plot on original age scale
#grid_age <- fhat$t * (max(swimmers$Age) - min(swimmers$Age)) + min(swimmers$Age)
#plot(grid_age, fhat$x, type = "l")

# mean estimation by pooling + splines
df_list <- df_norm |>
  group_by(ID) |>
  arrange(ID) |>
  group_split(.keep = TRUE) |>
  purrr::map(~list(t = .x$Age,
                   x = .x$Performance))

mu_spl <- splines_mean(X_list = df_list,
                       x = xout)

#plot(mu_spl$t, mu_spl$x, type = "l")

# estimation of eigenfunctions by zw
pca_ll <- fdapace::FPCA(Lt = purrr::map(df_list, ~.x$t),
                        Ly = purrr::map(df_list, ~.x$x),
                        optns = list(nRegGrid = length(xout),
                                     FVEthreshold = 0.95))

# plot(pca_ll$workGrid, pca_ll$phi[, 1], type = "l")
# plot(pca_ll$workGrid, pca_ll$phi[, 2], type = "l")
# plot(pca_ll$workGrid, pca_ll$phi[, 3], type = "l")

# Change the sign to negative for the 2nd and 3rd eigenfunction to make more
# sense - now they all look convex, but with different minimums
pca_ll$phi[, 2] <- -1 * pca_ll$phi[, 2]
pca_ll$phi[, 3] <- -1 * pca_ll$phi[, 3]


# interpolate estimated quantites onto observed design points
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
                       ~list(
                         t = df_list[[.x]]$t,
                         x = (df_list[[.x]]$x - mu_i[[.x]]) * phi_i[[.x]] / f_i[[.x]]
                         )
                       )

# cdf_hat <- list(t = xout,
#                 x = cdf_series(xi = df_age$Age,
#                                x = xout,
#                                J = 12,
#                                thresh = TRUE))

cdf_hat <- list(t = fhat$t,
                x = c(pracma::cumtrapz(x = fhat$t,
                                       y = fhat$x)))

# estimate scores using nn
xi_list <- purrr::map(varphi_i,
                      ~apply(.x$x, 2, function(varphi) {
                        mc_int(x = .x$t,
                               varphi = varphi,
                               cdf = cdf_hat)$varphi_int
                        })
                      )

df_mc <- do.call('rbind', xi_list)
colnames(df_mc) <- paste0("xi", seq_len(pca_ll$selectK))


# estimate scores using sample mean
xi_emp <- purrr::map(varphi_i, ~colMeans(.x$x))
df_emp <- do.call('rbind', xi_emp)
colnames(df_emp) <- paste0("xi", seq_len(pca_ll$selectK))

# estimate scores using riemann sums
xi_riemann <- purrr::map(varphi_i,
                         ~apply(.x$x, 2, function(xij) pracma::trapz(x = .x$t,
                                                                     y = xij)))

df_riemann <- do.call('rbind', xi_riemann)
colnames(df_riemann) <- paste0("xi", seq_len(pca_ll$selectK))



# scale to divide by empirical eigenvalues
df_scaled_mc <- scale(df_mc, center = FALSE)
df_scaled_emp <- scale(df_emp, center = FALSE)
df_scaled_riemann <- scale(df_riemann, center = FALSE)

# mc_kmeans <- kmeans(x = df_scaled_mc, centers = 2, nstart = 50)
# emp_kmeans <- kmeans(x = df_scaled_emp, centers = 2, nstart = 50)
# riemann_kmeans <- kmeans(x = df_scaled_riemann, centers = 2, nstart = 50)

mc_hclust <- hclust(dist(df_scaled_mc), method = "average")
cut_mc <- cutree(mc_hclust, k = 2)

emp_hclust <- hclust(dist(df_scaled_emp), method = "average")
cut_emp <- cutree(emp_hclust, k = 2)

riemann_hclust <- hclust(dist(df_scaled_riemann), method = "average")
cut_riemann <- cutree(riemann_hclust, k = 2)


# Analysis for clustering
id2_mc <- seq_along(cut_mc)[cut_mc == 2]
id2_emp <- seq_along(cut_emp)[cut_emp == 2]
id2_riemann <- seq_along(cut_riemann)[cut_riemann == 2]

ID_mc <- df_norm |>
  filter(ID %in% id2_mc) |>
  select(ID) |>
  distinct()

ID_emp <- df_norm |>
  filter(ID %in% id2_emp) |>
  select(ID) |>
  distinct()


ID_trapz <- df_norm |>
  filter(ID %in% id2_riemann) |>
  select(ID) |>
  distinct()

mc_emp_id_int <- intersect(ID_mc$ID, ID_emp$ID)
mc_trapz_id_int <- intersect(ID_mc$ID, ID_trapz$ID)

pal_manual_mc <- RColorBrewer::brewer.pal(n = nrow(ID_mc), name = "Paired")
names(pal_manual_mc) <- ID_mc$ID

pal_manual_emp <- RColorBrewer::brewer.pal(n = nrow(ID_emp), name = "Paired")
names(pal_manual_emp) <- ID_emp$ID=
pal_manual_emp[paste0(mc_emp_id_int)] <- pal_manual_mc[paste0(mc_emp_id_int)]

pal_manual_trapz <- RColorBrewer::brewer.pal(n = nrow(ID_trapz), name = "Paired")
names(pal_manual_trapz) <- ID_trapz$ID
pal_manual_trapz[paste0(mc_trapz_id_int)] <- pal_manual_mc[paste0(mc_trapz_id_int)]

L2_mc <- df_norm |>
  filter(ID %in% id2_mc) |>
  group_by(ID) |>
  mutate(ID = as.factor(ID)) |>
  ggplot(aes(x = Age, y = Performance, col = ID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  scale_color_manual(values = pal_manual_mc) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Individuals selected in cluster L2 (Control Neighbours)')

library(tikzDevice)
# tikz('cluster2_mc.tex',width=5,height=5, standAlone = TRUE)
# plot(L2_mc)
# dev.off()

# Set same colours for common curves

L2_emp <- df_norm |>
  filter(ID %in% id2_emp) |>
  group_by(ID) |>
  mutate(ID = as.factor(ID)) |>
  ggplot(aes(x = Age, y = Performance, col = ID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  scale_color_manual(values = pal_manual_emp) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Individuals selected in cluster L2 (Empirical Mean)')

# tikz('cluster2_emp.tex',width=5,height=5,standAlone = TRUE)
# plot(L2_emp)
# dev.off()

L2_riemann <- df_norm |>
  filter(ID %in% id2_riemann) |>
  group_by(ID) |>
  mutate(ID = as.factor(ID)) |>
  ggplot(aes(x = Age, y = Performance, col = ID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  scale_color_manual(values = pal_manual_trapz) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Individuals selected in cluster L2 (Trapezoidal Rule)')

# tikz('cluster2_riemann.tex',width=5,height=5,standAlone = TRUE)
# plot(L2_riemann)
# dev.off()



mu_range <- function(df, id) {
  df |>
    filter(ID %in% id) |>
    group_by(ID) |>
    summarise(range = max(Performance) - min(Performance)) |>
    ungroup() |>
    summarise(mu_range = mean(range))
}

mc_range <- sapply(list(seq_along(cut_mc)[cut_mc == 1], id2_mc), function(x) {
  mu_range(df = df_norm, id = x)$mu_range
})

riemann_range <- sapply(list(seq_along(cut_riemann)[cut_riemann == 1], id2_riemann), function(x) {
  mu_range(df = df_norm, id = x)$mu_range
})

emp_range <- sapply(list(seq_along(cut_riemann)[cut_emp == 1], id2_emp), function(x) {
  mu_range(df = df_norm, id = x)$mu_range
})


# plots

phi_df <- data.frame(phi1 = pca_ll$phi[, 1],
                     phi2 = pca_ll$phi[, 2],
                     phi3 = pca_ll$phi[, 3],
                     t = pca_ll$workGrid) |>
  pivot_longer(cols = phi1:phi3,
               names_to = "Eigenfunction",
               values_to = "phi_t")

phi_line <- phi_df |>
  ggplot(aes(x = t, y = phi_t, col = Eigenfunction))  +
  geom_line() +
  ylab("phi(t)") +
  ggtitle("Eigenfunctions of Performance Curves") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

tikz('phi_line.tex',width=5,height=5,standAlone = TRUE)
plot(phi_line)
dev.off()

tikz('mu_line.tex', width=5,height=5,standAlone = TRUE)
plot(mu_spl$t, mu_spl$x, type = "l",
     xlab = "t", ylab = "mu(t)", main = "Mean function (Performance Curves)")
dev.off()

tikz('f_line.tex', width=5,height=5,standAlone = TRUE)
plot(fhat$t, fhat$x, type = "l",
     xlab = "t", ylab = "f(t)", main = "Design Density (Performance Curves)")
dev.off()
