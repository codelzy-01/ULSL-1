#' Title
#'
#' @param Y the true label vector
#' @param predY the predicted label vector
#'
#' @return list containing ACC (accuracy), MIhat (normalized mutual information), and ARI
#' @export

ClusteringMeasure <- function(Y, predY) {

  if (ncol(as.matrix(Y)) != 1) {
    Y <- matrix(Y, ncol = 1)
  }
  if (ncol(as.matrix(predY)) != 1) {
    predY <- matrix(predY, ncol = 1)
  }

  n <- length(Y)

  # 重新编码标签从1开始连续编号
  uY <- unique(Y)
  nclass <- length(uY)
  Y0 <- numeric(n)
  if (nclass != max(Y)) {
    for (i in 1:nclass) {
      Y0[which(Y == uY[i])] <- i
    }
    Y <- Y0
  }

  uY <- unique(predY)
  nclass <- length(uY)
  predY0 <- numeric(n)
  if (nclass != max(predY)) {
    for (i in 1:nclass) {
      predY0[which(predY == uY[i])] <- i
    }
    predY <- predY0
  }

  # 计算纯度(Purity)
  Lidx <- unique(Y)
  classnum <- length(Lidx)
  predLidx <- unique(predY)
  pred_classnum <- length(predLidx)

  correnum <- 0
  for (ci in 1:pred_classnum) {
    incluster <- Y[which(predY == predLidx[ci])]
    if (length(incluster) > 0) {
      inclunub <- hist(incluster, breaks = seq(0.5, max(incluster)+0.5, 1), plot = FALSE)$counts
      correnum <- correnum + max(inclunub)
    }
  }
  Purity <- correnum / length(predY)

  # 使用最佳映射对齐标签
  res <- bestMap(Y, predY)

  # 计算准确率(ACC)

  ACC <- length(which(Y == res)) / length(Y)

  # 计算标准化互信息(MIhat)
  MIhat <- MutualInfo(Y, res)

  #计算ARI
  if (!requireNamespace("mclust", quietly = TRUE)) {
    install.packages("mclust")
  }
  ARI <- mclust::adjustedRandIndex(Y, res)

  return(list(ACC = ACC, MIhat = MIhat, Purity = Purity,ARI = ARI))
}
