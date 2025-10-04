#' Title
#'
#' @param H a d1 × n matrix (columns represent samples)
#' @param G  an n × n initial similarity matrix
#' @param beta regularization parameter
#' @param lambda feature difference penalty parameter
#'
#' @return an n × n similarity matrix Z with non-negative entries and rows summing to 1
#' @export

UpdateZ <- function(H, G, beta, lambda) {


  solveEta <- function(z_i) {
    sorted_z <- sort(z_i, decreasing = TRUE)
    s <- 0
    for (j in seq_along(sorted_z)) {
      s <- s + sorted_z[j]
      eta_prime <- (1 - s) / j
      if (eta_prime + sorted_z[j] <= 0) {
        eta_prime <- (1 - sum(sorted_z[1:(j - 1)])) / (j - 1)
        break
      }
    }
    if (!exists("eta_prime")) eta_prime <- (1 - sum(sorted_z)) / length(sorted_z)
    return(eta_prime)
  }

  n <- ncol(H)  # 样本数
  d1 <- nrow(H)


  HHt <- crossprod(H)  # H' * H，n x n
  H2 <- colSums(H^2)   # 每列的平方和
  D <- outer(H2, H2, "+") - 2 * HHt
  D <- lambda * D     # 加权

  Z <- matrix(0, nrow = n, ncol = n)

  for (i in 1:n) {
    z_i <- G[i, ] - D[i, ] / (2 * beta)
    eta_prime <- solveEta(z_i)
    Z[i, ] <- pmax(z_i + eta_prime, 0)
  }

  return(Z)
}
