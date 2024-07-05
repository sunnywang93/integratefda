library(mcscores)
library(parallel)
library(tictoc)
library(reshape2)
library(tidyr)
library(ggplot2)

n_basis <- 200
lambda_pow <- 2
x_length <- 51
ci_level <- 0.95
xi_cdf <- rt
#xi_cdf <- rchisq
#xi_cdf <- rnorm
points_dist <- runif
#cdf <- function(x) a/3*(x - 0.5)^3 + b*x + a/24
cdf <- identity
rout <- 500
a <- 3
b <- 1 - a/12
shift <- 0
sigma <- 0.1

x_cov <- seq(0, 1, length.out = 301)
# Generate large sample of curves to construct "true" covariance matrix
cov_pace <- purrr::map(seq_len(3000),
                     ~bm_kl_rd(k = n_basis, x = x_cov, lambda_rate = 2.8,
                               norm_constant = 0, norm_factor = 0)) |>
  purrr::map(~tcrossprod(.x$x)) |>
  (\(x) Reduce('+', x) / length(x))()


tic()
result <- mclapply(seq_len(rout), function(rep) {

x_obs <- sort(points_dist(n = x_length))
#x_obs <- sort(boat_cdf(a = a, n = x_length))
# x_list <- bm_kl_rd(k = n_basis, x = x_obs, lambda_rate = lambda_pow,
#                    norm_factor = sqrt(1/2), xi_dist = xi_cdf, df = 1)
#x_list <- kl_rd_poly(k = n_basis, x = x_obs, lambda_rate = lambda_pow)
x_list <- bm_kl_rd(k = n_basis, x = x_obs, lambda_rate = lambda_pow,
                   norm_constant = 0, norm_factor = sqrt(1/2), xi_dist = xi_cdf,
                   df = 4)

lag <- 1
eta <- rnorm(length(x_list$x) + lag)
em <- RcppRoll::roll_sum(x = eta, n = lag + 1) / sqrt(lag + 1)

#noise <- x_list$x * sigma * em
noise <- sigma * em
y_list <- list(t = x_obs,
               x = x_list$x,
               y = x_list$x + noise + shift)


lambda <- ((seq_len(n_basis) - 0.5) * pi)^(-lambda_pow)

mu_list <- list(t = x_obs,
                x = rep(0, length(x_obs)))

psi_list <- list(t = x_obs,
                 X = sqrt(2) * sin(outer(pi * x_obs, seq_len(n_basis) - 0.5)))


pace_noise <- 0.1


x_obs_bi <- expand.grid(y = x_obs,
                        x = x_obs)

cov_pace_rd <- pracma::interp2(x = x_cov,
                               y = x_cov,
                               Z = cov_pace,
                               xp = x_obs_bi[, 2],
                               yp = x_obs_bi[, 1]) |>
  matrix(nrow = length(x_obs),
         ncol = length(x_obs))

Sigma_mat <- cov_pace_rd + pace_noise * diag(nrow = nrow(cov_pace_rd),
                                             ncol = ncol(cov_pace_rd))


pdf_list <- list(t = x_obs,
                 x = rep(1, length(x_obs)))

xi_hat <- mc_scores(y_list = y_list,
                    mu_list = mu_list,
                    psi_list = psi_list,
                    pdf_list = pdf_list,
                    cdf = cdf,
                    noise_sd = sigma)

ci_exact <- confint_mc(xi_hat, ci_level)
ci_lim <- confint_lim(xi_hat, pdf_list, ci_level)

xi_pace <- pace_scores(y_list, lambda, psi_list, Sigma_mat, mu_list)

ci_pace <- confint_pace(scores_vec = xi_pace,
                        lambda_vec = lambda,
                        phi_list = psi_list,
                        Sigma_mat = Sigma_mat,
                        conf_level = ci_level)

if(shift != 0) {
  xi_true <- x_list$xi + shift * sqrt(2) / ((seq_len(n_basis) - 0.5) * pi)
} else {
  xi_true <- x_list$xi
}



exact_count <- (xi_true >= ci_exact$ci_l) & (xi_true <= ci_exact$ci_u)
limit_count <- (xi_true >= ci_lim$ci_l) & (xi_true <= ci_lim$ci_u)
pace_count <- (xi_true >= ci_pace$ci_l) & (xi_true <= ci_pace$ci_u)

xi_int <- apply(psi_list$X, 2, function(x) pracma::trapz(x = psi_list$t,
                                                         y = x * x_list$x))

list(exact_count = exact_count,
     lim_count = limit_count,
     pace_count = pace_count,
     width_exact = ci_exact$ci_u - ci_exact$ci_l,
     width_lim = ci_lim$ci_u - ci_lim$ci_l,
     width_pace = ci_pace$ci_u - ci_pace$ci_l,
     mc_err = abs(xi_hat$xi_hat - xi_true),
     bias_mc = xi_hat$xi_hat - xi_true,
     bias_pace = xi_pace - xi_true,
     pace_err = abs(xi_pace - xi_true),
     xi_mc = xi_hat$xi_hat,
     xi_pace = xi_pace,
     xi_true = xi_true,
     ci_exact = ci_exact,
     ci_pace = ci_pace)

}, mc.cores = 5)
toc()


# Mean Absolute Errors =========================================================
mc_err_mat <- t(sapply(result, function(r) r$mc_err[1:5]))
colnames(mc_err_mat) <- paste0("xi", seq_len(5))
pace_err_mat <- t(sapply(result, function(r) r$pace_err[1:5]))
colnames(pace_err_mat) <- paste0("xi", seq_len(5))

df1 <- as.data.frame(mc_err_mat)
df2 <- as.data.frame(pace_err_mat)

df1$Source <- "mc"
df2$Source <- "pace"

combined_df <- rbind(df1, df2)

melted_df <- melt(combined_df, id.vars = "Source")

ggplot(melted_df, aes(x = variable, y = value, fill = Source)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots for scores", x = "Scores", y = "MAE")




# Check bias vs variance ======================================================
bias_mc <- t(sapply(result, function(r) r$bias_mc[1:5]))
colnames(bias_mc) <- paste0("xi", seq_len(5))
bias_pace <- t(sapply(result, function(r) r$bias_pace[1:5]))
colnames(bias_pace) <- paste0("xi", seq_len(5))

df1 <- as.data.frame(bias_mc)
df2 <- as.data.frame(bias_pace)

# Add an identifier column to each data frame
df1$Source <- "mc"
df2$Source <- "pace"

# Combine the two data frames
combined_df <- rbind(df1, df2)

# Melt the combined data frame for ggplot2
melted_df <- melt(combined_df, id.vars = "Source")

# Create the boxplot
ggplot(melted_df, aes(x = variable, y = value, fill = Source)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots for scores", x = "Scores", y = "MAE")



# Check coverage ================================================================
exact_coverage <- Reduce('+', purrr::map(result, ~.x$exact_count))
lim_coverage <- Reduce('+', purrr::map(result, ~.x$lim_count))
pace_coverage <- Reduce('+', purrr::map(result, ~.x$pace_count))



