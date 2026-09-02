#' Euclidean distance from each observation to a grid point
#' @keywords internal
#' @noRd
euc.dist <- function(biom, x) {
  xmat <- do.call("rbind", replicate(nrow(biom), x, simplify = FALSE))
  sqrt(rowSums((biom - xmat)^2))
}
