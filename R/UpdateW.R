#' @title updated matrix W
#'
#' @param Coef latent representation matrix H
#' @param Data data matrix M
#' @param D_Mat the iteratively updated W matrix
#'
#' @return the updated matrix W
#' @export

UpdateW <- function(Coef, Data, D_Mat) {
  # Coef 表示 H，Data 表示 M，D_Mat 表示之前的 W

  Imat <- diag(nrow(Coef))  # 单位矩阵，行数与 Coef 相同

  TempCoef <- as.matrix(Coef)
  TempData <- as.matrix(Data)
  rho <- 1
  rate_rho <- 1.2
  TempS <- D_Mat
  TempT <- matrix(0, nrow = nrow(TempS), ncol = ncol(TempS))
  previousD <- D_Mat
  Iter <- 1
  ERROR <- 1

  while (ERROR > 1e-6 && Iter < 100) {
    # 计算 TempD
    TempD <- (rho * (TempS - TempT) + TempData %*% t(Coef)) %*% solve(rho * Imat + Coef %*% t(Coef))
    # 计算 TempS，进行列归一化，并确保列范数 <= 1
    TempS <- normcol_lessequal(TempD+TempT)

    # 更新 TempT
    TempT <- TempT + TempD - TempS

    # 更新 rho
    rho <- rate_rho * rho

    # 误差计算
    ERROR <- mean((previousD - TempD)^2)

    # 更新 previousD
    previousD <- TempD

    Iter <- Iter + 1
  }

  return(TempD)
}
