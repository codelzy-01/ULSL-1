#' @title Unified Latent and Similarity Learning  (ULSL)
#' @description
#' This function implements the ULSL algorithm, which integrates multiple
#' omics views (or data matrices) into a unified latent representation
#' and similarity graph, and performs clustering.
#' It uses an iterative optimization approach to jointly learn:
#' \itemize{
#'   \item Projection matrix \( W \)
#'   \item Latent representation matrix \( H \)
#'   \item Similarity matrix \( Z \)
#' }
#' The final output includes cluster assignments and optionally clustering performance measures.
#'
#' @param X A list of numeric matrices, each representing a view of the data
#'          (features × samples).
#'          Each matrix is expected to have the same number of columns (samples).
#' @param sss Numeric, regularization parameter controlling trade-off in objective function (default: 0.5).
#' @param alpha0 Numeric, base regularization parameter for similarity alignment (default: 2).
#' @param beta Numeric, regularization parameter for Laplacian constraint (default: 0.1).
#' @param d Integer, latent dimension size (default: 40).
#' @param maxIters Integer, maximum number of optimization iterations (default: 200).
#' @param gt Optional numeric or factor vector, ground-truth cluster labels.
#'           If provided, the function returns clustering performance measures.
#' @param c Optional integer, number of clusters.
#'          If not provided, the algorithm automatically determines it.
#' @param custom_range Optional numeric vector of length 2 specifying candidate cluster range.
#'                     If NULL, it will be determined automatically based on sample size.
#'
#' @return
#' If `gt` is provided:
#' \item{measures}{A list containing clustering evaluation metrics (Accuracy, NMI, Purity, etc.).}
#' Otherwise:
#' \item{cl}{A vector of cluster assignments.}
#' \item{Obj}{A numeric vector of objective function values across iterations.}
#' \item{Z}{The final learned similarity matrix.}
#'
#' @details
#' ULSL proceeds in the following steps:
#' \enumerate{
#'   \item Normalize each view to ensure comparability.
#'   \item Compute pairwise distances and affinity matrices for each view.
#'   \item Fuse affinity matrices using the SNF (Similarity Network Fusion) algorithm.
#'   \item Initialize latent representations \( W \) and \( H \).
#'   \item Iteratively update \( W \), \( H \), and \( Z \) until convergence or reaching `maxIters`.
#'   \item If the number of clusters `c` is not given:
#'         \itemize{
#'           \item Determine candidate cluster numbers (`NUMC`) from `custom_range` or automatically.
#'           \item Select the optimal number of clusters using `nemo.num.clusters`.
#'         }
#'   \item Perform spectral clustering on \( Z \) to obtain cluster labels.
#'   \item If `gt` is provided, compute clustering evaluation metrics.
#' }
#'

#' @export

ULSL <- function(X, sss =0.5 , alpha0 =2 , beta =0.1 , d =40 , maxIters =200 , gt = NULL, c = NULL,custom_range = NULL) {


  if (!is.null(gt)) {
    c <- length(unique(gt))
  } else {
    if (is.null(c)) {
      #stop("You must provide either 'gt' or a value for 'c'.")
    }
  }
  library(SNFtool)

  T <- 30
  V <- length(X)
  N <- ncol(X[[1]])  # 数据点数量
  K <- N/10
  MAXiter <- 1000
  REPlic <- 200
  #omega <- rep(1 / V, V)

  zr <- 1e-10

  D1 <- sapply(X, function(x) nrow(x))  # 每个视图的特征维度

  X_guiyi <- list()
  for (i in 1:V) {
    X_guiyi[[i]] <- X[[i]] + .Machine$double.eps
    row_sums <- sqrt(rowSums(X_guiyi[[i]]^2))
    X_guiyi[[i]] <- as.matrix(sweep(X_guiyi[[i]], 1, row_sums, "/"))
  }






  Dist <- vector("list", V)
  for (i in 1:V) {
    Dist[[i]] <- (dist2(t(as.matrix(X_guiyi[[i]])),t(as.matrix(X_guiyi[[i]]))))
  }


  S <- lapply(Dist, function(everyd) affinityMatrix(everyd, K, 0.5))


  A <- SNF(S, K, T)
  A <- (A + t(A)) / 2
  DD <- diag(rowSums(A))
  L <- DD - A

  M <- do.call(rbind, X_guiyi)
  M <- as.matrix(M)

  alpha <- alpha0*(dim(M)[1]/N)


  SD <- sum(D1)
  W <- matrix(0, nrow = SD, ncol = d)

  #set.seed(123)
  H <- matrix(runif(d * N), nrow = d, ncol = N)



  Obj <- numeric(maxIters)

  for (it in 1:maxIters) {
    # 更新 W（需要定义 UpdateW 函数）
    W <- UpdateW(H, M, W)

    # 更新 H（需要定义 UpdateH 函数）
    H <- UpdateH(M, W, L, beta)

    #更新Z
    Z <- UpdateZ(H, A, alpha, beta)
    Z_new <- Z
    # 更新 L（基于当前 Z）
    Z <- (Z + t(Z)) / 2
    DZ <- diag(rowSums(Z))
    L <- DZ - Z

    # 目标函数
    Obj[it] <- norm(M - W %*% H, "F")^2 + beta * sum(diag(H %*% L %*% t(H))) + alpha *norm(Z_new - A, "F")^2
    if (it > 1 && abs(Obj[it] - Obj[it - 1]) / Obj[it - 1] < 1e-4) {
      break
    }
  }


  print(norm(M - W %*% H, "F")^2)
  print(beta * sum(diag(H %*% L %*% t(H))))
  print(alpha *norm(Z_new - A, "F")^2)
  # 自动确定聚类数 C


  if (is.null(c) && is.null(gt)) {

    if (is.null(custom_range)) {
      NUMC <- get_cluster_range(N, custom_range)
    } else {
      NUMC <- custom_range
    }

    c <- num.clusters(Z, NUMC)
    cl <- spectralClustering(Z, c)

  } else {
    cl <- spectralClustering(Z,c)
  }

  # 聚类

  cl
  table(cl)
  # 最后的部分逻辑判断
  if (!is.null(gt)) {
    measures <- ClusteringMeasure(gt, cl)
    return(measures)
  } else {
    return(list(cl = cl, Obj = Obj,Z=Z))
  }


}
