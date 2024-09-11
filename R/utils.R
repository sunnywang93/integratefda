#' Compute the scores using conditional expectation of PACE
#'
#' @param y_list List of the observations, containing the elements:
#' - **$t** Vector of sampling points.
#' - **$y** Vector of observed points.
#' @param lambda_vec Vector containing the eigenvalues.
#' @param phi_list List of the eigenfunctions, containing the elements:
#' - **$t** Vector of sampling points.
#' - **$X** Matrix of observed points, with sampling points on rows and the index
#' @param Sigma_mat Matrix containing the covariance plus noise.
#' @param mu_list List of the mean function, containing the elements:
#' - **$t** Vector of sampling points.
#' - **$x** Vector of observed points.
#' @returns Vector containing the estimated scores.
#' @export


pace_scores <- function(y_list, lambda_vec, phi_list, Sigma_mat, mu_list) {

  yi_cen <- y_list$y - mu_list$x

  phi_norm <- sweep(phi_list$X, 2, lambda_vec, FUN = "*")

  sapply(seq_len(length(lambda_vec)),
         function(j) crossprod(phi_norm[, j], solve(Sigma_mat)) %*% yi_cen)

}

#' Compute confidence intervals based on the PACE method
#'
#' @param scores_vec Vector containing the scores.
#' @param lambda_vec Vector containing the eigenvalues.
#' @param phi_list List of the eigenfunctions, containing the elements:
#' - **$t** Vector of sampling points.
#' - **$X** Matrix of observed points, with sampling points on rows and the index.
#' @param Sigma_mat Matrix containing the covariance plus noise.
#' @param conf_level Numeric, indicating the confidence level in which the
#' intervals should be constructed.
#' @returns List, containing the elements:
#' - **$scores** Vector, containing the estimated scores.
#' - **$ci_l** Vector, containing the lower bound of the confidence interval.
#' - **$ci_u** Vector, containing the upper bound of the confidence interval.
#' - **$sm** Vector, containing the standard deviation.
#' @export

confint_pace <- function(scores_vec, lambda_vec, phi_list, Sigma_mat, conf_level) {

  H_mat <- t(sweep(phi_list$X, 2, lambda_vec, FUN = "*"))

  omega_k <- diag(lambda_vec) - crossprod(t(H_mat),
                                          tcrossprod(solve(Sigma_mat), H_mat))

  sm <- sqrt(diag(omega_k))

  crit <- 1 - conf_level

  ci_l <- scores_vec - qnorm(1 - crit/2) * sm
  ci_u <- scores_vec + qnorm(1 - crit/2) * sm

  list(scores = scores_vec,
       ci_l = ci_l,
       ci_u = ci_u,
       sm = sm)

}

#' Computes metric for binary classifiers
#'
#' @param y_hat Vector, the predicted classes.
#' @param y_true Vector, the true classes.
#' @param type String, either `f1` or `accuracy`.
#' @returns Numeric.
#' @export

class_score <- function(y_hat, y_true, type = "f1") {

  conf_matrix <- table(y_hat, y_true)
  tp <- conf_matrix[2, 2]
  fp <- conf_matrix[2, 1]
  fn <- conf_matrix[1, 2]
  tn <- conf_matrix[1, 1]

  precis <- tp / (tp + fp)
  recall <- tp / (tp + fn)

  if(type == "f1") {
    2 * precis * recall / (precis + recall)
  } else {
    (tp + tn) / (tp + tn + fp + fn)
  }

}




