#' Title
#'
#' @param M data matrix M
#' @param W the iteratively updated W matrix
#' @param L Laplacian matrix
#' @param lambda feature difference penalty parameter
#'
#' @return the updated matrix H
#' @export

UpdateH <- function(M, W, L, lambda) {
  # 加载 numpy 和 scipy.linalg
  np <- reticulate::import("numpy")
  scipy <- reticulate::import("scipy.linalg")

  # 将 R 矩阵转换为 numpy 数组
  M_np <- np$array(M)
  W_np <- np$array(W)
  L_np <- np$array(L)

  # 获取样本数量 N
  N <- dim(M)[2]

  # 计算矩阵 E = W'W 和 Q = W'M
  E <- np$dot(np$transpose(W_np), W_np)
  Q <- np$dot(np$transpose(W_np), M_np)

  # 定义 A, B, C
  A <- np$transpose(E)
  B <- lambda * np$transpose(L_np)
  C <- Q

  # 使用 solve_sylvester 求解 Sylvester 方程: A·H + H·B = C
  H_np <- scipy$solve_sylvester(A, B, C)

  # 将 Python 的结果转换为 R 矩阵返回
  H <- reticulate::py_to_r(H_np)
  return(H)
}
