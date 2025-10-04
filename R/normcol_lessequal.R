#' @title Normalize Columns of a Matrix to Unit Length (If Necessary)
#' @description
#' This function normalizes each column of the input matrix so that
#' its Euclidean norm does not exceed 1.
#' If a column already has a norm less than or equal to 1,
#' it remains unchanged.
#'
#' @param matin A numeric matrix whose columns are to be normalized.
#'
#' @return A numeric matrix with columns normalized such that
#'         their Euclidean norms are ≤ 1.
#'
#'
#' @export

normcol_lessequal <- function(matin) {
  col_norms <- sqrt(colSums(matin^2)) + .Machine$double.eps
  matin / rep(pmax(1, col_norms), each = nrow(matin))
}
