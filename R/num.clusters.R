#' @title Estimate the Optimal Number of Clusters via Eigengap Analysis
#' @description
#' This function estimates the optimal number of clusters in a dataset
#' based on the eigengap method from spectral clustering.
#' It computes the normalized graph Laplacian from the given similarity matrix,
#' calculates eigenvalues, and selects the number of clusters that maximizes
#' the scaled eigengap within a candidate range.
#'
#' @param W A symmetric similarity matrix (affinity matrix) representing relationships between samples.
#' @param NUMC An integer vector of candidate cluster numbers to evaluate (default: 2:10).
#'
#' @return An integer representing the estimated optimal number of clusters.

#' @export
num.clusters <- function(W, NUMC=2:10) {
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC = NUMC[NUMC > 1]
  }
  W = (W + t(W))/2
  diag(W) = 0
  if (length(NUMC) > 0) {
    degs = rowSums(W)
    degs[degs == 0] = .Machine$double.eps
    D = diag(degs)
    L = D - W
    Di = diag(1/sqrt(degs))
    L = Di %*% L %*% Di
    print(dim(L))
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return = T)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]
    eigengap = abs(diff(eigs$values))
    eigengap = (1:length(eigengap)) * eigengap




    t1 <- sort(eigengap[NUMC], decreasing = TRUE, index.return = T)$ix
    t1
    return(NUMC[t1[1]])
  }
}
