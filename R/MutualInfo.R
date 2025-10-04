#' Title
#'
#' @param L1 the first label vector
#' @param L2 the second label vector
#'
#' @return the normalized mutual information (NMI) value
#' @export
#'

MutualInfo <- function(L1, L2) {


  L1 <- matrix(L1, ncol = 1)
  L2 <- matrix(L2, ncol = 1)

  if (length(L1) != length(L2)) {
    stop("size(L1) must == size(L2)")
  }

  # 将标签调整为从1开始
  L1 <- L1 - min(L1) + 1
  L2 <- L2 - min(L2) + 1

  # 创建联合分布矩阵
  nClass <- max(max(L1), max(L2))
  G <- matrix(0, nrow = nClass, ncol = nClass)

  for (i in 1:nClass) {
    for (j in 1:nClass) {
      G[i, j] <- length(which(L1 == i & L2 == j)) + .Machine$double.eps
    }
  }

  sumG <- sum(G)

  # 计算边际分布
  P1 <- rowSums(G) / sumG
  P2 <- colSums(G) / sumG

  # 计算熵
  H1 <- -sum(P1 * log2(P1))
  H2 <- -sum(P2 * log2(P2))

  # 计算联合分布和互信息
  P12 <- G / sumG
  PPP <- P12 / outer(P1, P2)
  PPP[abs(PPP) < 1e-12] <- 1
  MI <- sum(P12 * log2(PPP))

  # 标准化互信息
  MIhat <- MI / max(H1, H2)

  return(Re(MIhat))  # 确保返回实数
}
