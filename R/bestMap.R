#' Title
#'
#' @param L1 the true label vector
#' @param L2 the predicted label vector
#'
#' @return label vector
#' @export
#' @description
#' Use the Hungarian algorithm to optimally adjust the L2 label assignments so that they best match those of L1.

bestMap <- function(L1, L2) {


  if (length(L1) != length(L2)) {
    stop("L1 和 L2 的长度必须一致")
  }

  # 将标签标准化为连续正整数
  L1 <- as.integer(factor(L1))
  L2 <- as.integer(factor(L2))

  nClass <- max(max(L1), max(L2))

  # 构造混淆矩阵
  G <- table(L1, L2)
  G <- as.matrix(G)

  # 使用匈牙利算法（注意转置）
  if (!requireNamespace("clue", quietly = TRUE)) {
    install.packages("clue")
  }
  library(clue)

  matching <- solve_LSAP(t(G), maximum = TRUE)  # 转置！列→行匹配
  newL2 <- matching[L2]  # 重新映射 L2 标签

  return(as.integer(newL2))
}
