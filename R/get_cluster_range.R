#' @title Determine Range of Candidate Cluster Numbers
#' @description
#' This function determines an appropriate range of candidate cluster numbers
#' for clustering analysis, based either on a user-specified range or
#' automatically depending on the sample size.
#'
#' @param sample_size An integer specifying the number of samples in the dataset.
#' @param custom_range An optional numeric vector of length 2 specifying
#'                     the lower and upper bounds of the cluster range.
#'                     If provided, this range is returned directly.
#'
#' @return An integer sequence representing the range of candidate cluster numbers.
#'
#' @export
get_cluster_range <- function(sample_size, custom_range = NULL) {

  if (!is.null(custom_range)) {
    stopifnot(length(custom_range) == 2, custom_range[1] <= custom_range[2])
    return(seq(custom_range[1], custom_range[2]))
  }

  # 根据样本量设置默认范围
  if (sample_size > 500) {
    return(4:8)
  } else if (sample_size >= 200) {
    return(3:6)
  } else {
    return(2:5)
  }
}
