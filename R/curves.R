#' Recover the curves using the KL-decomposition
#'
#' @param scores_vec Vector of scores
#' @param psi_list List of the eigenfunctions, containing the elements:
#' - **$t** Vector of sampling points.
#' - **$X** Matrix of observed points, with sampling points on rows and the index
#' of the eigenfunctions on the columns.
#' @returns List containing the sampling points and observed points.


kl_curve <- function(scores_vec, psi_list) {

  X <- sweep(psi_list$X, 2, scores_vec, FUN = "*") |>
    rowSums()

  list(t = psi_list$t,
       x = X)

}














