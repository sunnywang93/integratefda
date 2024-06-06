sigma <- 0.1
n_basis <- 200

test <- lapply(seq_len(500), function(rep) {


t_unif <- sort(runif(n = 101))
x_list <- bm_kl_rd(k = n_basis, x = t_unif)
y_list <- list(t = t_unif,
               x = x_list$x,
               y = x_list$x + rnorm(length(t_unif), 0, sigma))


mu_list <- list(t = t_unif,
                x = rep(0, length(t_unif)))


psi_list <- list(t = t_unif,
                 X = sqrt(2) * sin(outer(pi * t_unif, seq_len(n_basis) - 0.5)))

pdf_list <- list(t = t_unif,
                 x = rep(1, length(t_unif)))

cdf <- identity

xi_hat <- mc_scores(y_list = y_list,
                    mu_list = mu_list,
                    psi_list = psi_list,
                    pdf_list = pdf_list,
                    cdf = cdf,
                    noise_sd = sigma)

mean(abs(xi_hat$xi_hat - x_list$xi))

ci_exact <- confint_mc(xi_hat, 0.95)
ci_lim <- confint_lim(xi_hat, pdf_list, 0.95)

par(mfrow = c(1, 2))
qqnorm((xi_hat$xi_hat - x_list$xi) / ci_exact$sm, main = "Cond Var")
qqline((xi_hat$xi_hat - x_list$xi) / ci_exact$sm, col = "steelblue", pch = 2)

qqnorm((xi_hat$xi_hat - x_list$xi) / ci_lim$sm, main = "Lim Var")
qqline((xi_hat$xi_hat - x_list$xi) / ci_lim$sm, col = "steelblue", pch = 2)


exact_count <- sapply(seq_along(x_list$xi),
       function(id) (x_list$xi[id] >= ci_exact$ci_l[id]) &
         (x_list$xi[id] <= ci_exact$ci_u[id]))

limit_count <- sapply(seq_along(x_list$xi),
       function(id) (x_list$xi[id] >= ci_lim$ci_l[id]) &
         (x_list$xi[id] <= ci_lim$ci_u[id]))


list(exact_count = exact_count,
     lim_count = limit_count)


})

exact_pct <- Reduce('+', purrr::map(test, ~.x$exact_count))
lim_pct <- Reduce('+', purrr::map(test, ~.x$lim_count))

exact_pct[1:10] / 500
exact_pct[1:10] / 500
