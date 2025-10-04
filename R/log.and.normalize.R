#' Log-transform and normalize multi-omics data
#'
#' This function performs optional log-transformation and normalization
#' on a list of omics datasets. Each dataset is assumed to be a matrix
#' with rows representing features (e.g., genes) and columns representing samples.
#'
#' The function includes three main preprocessing steps:
#' \enumerate{
#'   \item Remove features with zero variance.
#'   \item Apply \code{log(1+x)} transformation (optional).
#'   \item Normalize each dataset (optional).
#' }
#'
#' @param omics.data A list of omics datasets. Each element should be a numeric
#'   matrix with features in rows and samples in columns.
#' @param normalize Logical or logical vector. If a single value, it is recycled
#'   to all views. If a vector, its length must match the number of datasets.
#'   Indicates whether to apply normalization to each dataset.
#' @param log.transform Logical or logical vector. If a single value, it is
#'   recycled to all views. If a vector, its length must match the number of
#'   datasets. Indicates whether to apply \code{log(1+x)} transformation
#'   to each dataset.
#'
#' @return A list of processed omics datasets with the same structure as
#'   the input \code{omics.data}.
#'
#' @export
log.and.normalize <- function(omics.data,
                              normalize = TRUE,
                              log.transform = FALSE) {
  n_views <- length(omics.data)

  # 如果 normalize 或 log.transform 是单一值，扩展为长度为视图数量的向量
  if (length(normalize) == 1) normalize <- rep(normalize, n_views)
  if (length(log.transform) == 1) log.transform <- rep(log.transform, n_views)

  # 检查长度
  if (length(normalize) != n_views) {
    stop("normalize must be either length 1 or equal to the number of views (omics.data)")
  }
  if (length(log.transform) != n_views) {
    stop("log.transform must be either length 1 or equal to the number of views (omics.data)")
  }

  # 1. 过滤掉无方差特征
  for (i in seq_along(omics.data)) {
    omics.data[[i]] <- omics.data[[i]][apply(omics.data[[i]], 1, var) > 0, ]
  }

  # 2. 按需做 log(1 + x) 变换
  for (i in seq_along(omics.data)) {
    if (log.transform[i]) {
      omics.data[[i]] <- log(1 + omics.data[[i]])
    }
  }

  # 3. 归一化
  for (i in seq_along(omics.data)) {
    if (normalize[i]) {
      omics.data[[i]] <- normalize.matrix(omics.data[[i]])
    }
  }

  return(omics.data)
}
