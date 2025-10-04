#' @title Hungarian Algorithm for Optimal Assignment
#' @description
#' Implements a variant of the Hungarian algorithm to find an optimal assignment
#' in a cost matrix. This is often used for problems such as label matching,
#' where one needs to align two sets of labels with minimal cost.
#' The function returns the updated cost matrix, the assignment vector,
#' and a list of unassigned rows.
#'
#' @param A A numeric matrix representing the cost matrix.
#'          Negative values indicate potential matches between rows and columns.
#'
#' @return A list containing:
#' \describe{
#'   \item{A}{The updated cost matrix after performing the assignment process.}
#'   \item{C}{An integer vector of length equal to the number of columns in A,
#'            indicating the row index assigned to each column. A value of 0
#'            means the column is unassigned.}
#'   \item{U}{An integer vector containing the indices of unassigned rows,
#'            with the first element reserved as a placeholder.}
#' }

#'

#'
#' @export

hminiass <- function(A) {
  n <- nrow(A)
  np1 <- ncol(A)

  # 初始化返回向量
  C <- numeric(n)
  U <- numeric(n+1)

  # 初始化零元素指针
  LZ <- numeric(n)
  NZ <- numeric(n)

  for (i in 1:n) {
    # 获取行i中第一个未分配的零
    lj <- np1
    j <- -A[i, lj]

    # 循环直到找到未分配列或列表结束
    while (j != 0 && C[j] != 0) {
      lj <- j
      j <- -A[i, lj]
    }

    if (j != 0) {
      # 找到未分配列的零
      C[j] <- i
      A[i, lj] <- A[i, j]
      NZ[i] <- -A[i, j]
      LZ[i] <- lj
      A[i, j] <- 0
    } else {
      # 检查当前行所有零
      lj <- np1
      j <- -A[i, lj]

      while (j != 0) {
        r <- C[j]
        lm <- LZ[r]
        m <- NZ[r]

        while (m != 0) {
          if (C[m] == 0) break
          lm <- m
          m <- -A[r, lm]
        }

        if (m == 0) {
          lj <- j
          j <- -A[i, lj]
        } else {
          # 替换零元素
          A[r, lm] <- -j
          A[r, j] <- A[r, m]
          NZ[r] <- -A[r, m]
          LZ[r] <- j
          A[r, m] <- 0
          C[m] <- r

          # 移除A[i,j]
          A[i, lj] <- A[i, j]
          NZ[i] <- -A[i, j]
          LZ[i] <- lj
          A[i, j] <- 0
          C[j] <- i
          break
        }
      }
    }
  }

  # 创建未分配行列表
  empty <- which(C == 0)
  U <- numeric(n+1)
  U[c(n+1, empty)] <- c(empty, 0)

  return(list(A = A, C = C, U = U))
}
