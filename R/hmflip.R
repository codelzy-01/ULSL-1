#' @title Path-Flipping Step in the Hungarian Algorithm
#' @description
#' This function implements the path-flipping step used in a variant of
#' the Hungarian algorithm for optimal assignment problems.
#' Given a cost matrix and current assignment state, it updates the
#' assignment along a path to improve the overall matching.
#'
#' @param A A numeric matrix representing the current cost matrix,
#'          where negative values are used as pointers for unassigned zeros.
#' @param C An integer vector indicating the current column-to-row assignment.
#' @param LC An integer vector storing the assigned row for each column.
#' @param LR An integer vector storing the assigned column for each row.
#' @param U An integer vector of unassigned rows, with U reserved as a placeholder.
#' @param l An integer representing the current column index in the flipping path.
#' @param r An integer representing the current row index in the flipping path.
#'
#' @return A list containing:
#' \describe{
#'   \item{A}{The updated cost matrix after the flipping step.}
#'   \item{C}{The updated column-to-row assignment vector.}
#'   \item{U}{The updated unassigned row list.}
#' }
#'

#'

#'
#' @export
hmflip <- function(A, C, LC, LR, U, l, r) {
  n <- nrow(A)

  while (TRUE) {
    # 在列l中分配行r
    C[l] <- r

    # 从零列表中移除该零
    m <- which(A[r,] == -l)
    A[r, m] <- A[r, l]
    A[r, l] <- 0

    # 如果是路径的第一个零
    if (LR[r] < 0) {
      U[n+1] <- U[r]
      U[r] <- 0
      return(list(A = A, C = C, U = U))
    } else {
      # 回溯路径
      l <- LR[r]
      A[r, l] <- A[r, ncol(A)]
      A[r, ncol(A)] <- -l
      r <- LC[l]
    }
  }
}
