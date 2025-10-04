#' @title Reduction Step in the Hungarian Algorithm
#' @description
#' This function implements the reduction step in a variant of the Hungarian algorithm
#' for optimal assignment problems. It adjusts the cost matrix by reducing uncovered
#' elements and adjusting covered elements, thereby creating new zeros
#' and enabling further assignments.
#'
#' @param A A numeric matrix representing the current cost matrix,
#'          where negative values indicate pointers for unassigned zeros.
#' @param CH An integer vector representing the last assigned zero for each row.
#' @param RH An integer vector storing rows that need to be explored.
#' @param LC An integer vector representing the current column assignments.
#' @param LR An integer vector representing the current row assignments.
#' @param SLC An integer vector of columns to be scanned (covered columns list).
#' @param SLR An integer vector of rows to be scanned (covered rows list).
#'
#' @return A list containing:
#' \describe{
#'   \item{A}{The updated cost matrix after the reduction step.}
#'   \item{CH}{The updated list of last assigned zeros for each row.}
#'   \item{RH}{The updated list of rows to be explored.}
#' }

#' @export

hmreduce <- function(A, CH, RH, LC, LR, SLC, SLR) {
  n <- nrow(A)

  # 找出未覆盖的行和列
  coveredRows <- LR == 0
  coveredCols <- LC != 0
  r <- which(!coveredRows)
  c <- which(!coveredCols)

  # 找出未覆盖元素的最小值
  m <- min(A[r, c])

  # 从未覆盖元素中减去最小值
  A[r, c] <- A[r, c] - m

  # 检查所有未覆盖的列
  for (j in c) {
    # 按路径顺序检查未覆盖的行
    for (i in SLR) {
      # 如果是新的零
      if (A[i, j] == 0) {
        # 如果行不在未探索列表中
        if (RH[i] == 0) {
          RH[i] <- RH[n+1]
          RH[n+1] <- i
          CH[i] <- j
        }

        # 找出行i中最后一个未分配的零
        row <- A[i,]
        colsInList <- -row[row < 0]
        if (length(colsInList) == 0) {
          l <- ncol(A)
        } else {
          l <- colsInList[row[colsInList] == 0]
          if (length(l) == 0) l <- ncol(A)
        }

        # 将该零添加到列表末尾
        A[i, l] <- -j
      }
    }
  }

  # 给双重覆盖的元素加上最小值
  r <- which(coveredRows)
  c <- which(coveredCols)
  A[r, c] <- A[r, c] + m

  # 处理要移除的零
  zeros <- which(A[r, c] <= 0, arr.ind = TRUE)
  if (length(zeros) > 0) {
    i <- r[zeros[,1]]
    j <- c[zeros[,2]]

    for (k in 1:length(i)) {
      lj <- which(A[i[k],] == -j[k])
      A[i[k], lj] <- A[i[k], j[k]]
      A[i[k], j[k]] <- 0
    }
  }

  return(list(A = A, CH = CH, RH = RH))
}
