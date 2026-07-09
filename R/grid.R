#' Build a biomarker evaluation grid
#'
#' Creates an evenly spaced grid over the observed range of each biomarker
#' column, using `expand.grid()` for multi-dimensional combinations.
#'
#' @param biom Biomarker values: a vector, matrix, or data.frame.
#' @param grid.size Number of grid points per biomarker dimension.
#' @return A data.frame of grid coordinates with the same column names as `biom`.
#' @export
#' @examples
#' grid_biom(data.frame(x = c(1, 5), y = c(2, 8)), grid.size = 3)
grid_biom <- function(biom, grid.size = 25) {
  dbiom <- as.data.frame(biom)
  grid.cut <- lapply(seq_len(ncol(dbiom)), function(j) {
    seq(min(dbiom[[j]], na.rm = TRUE), max(dbiom[[j]], na.rm = TRUE),
        length.out = grid.size)
  })
  names(grid.cut) <- names(dbiom)
  as.data.frame(expand.grid(grid.cut))
}

#' @keywords internal
#' @noRd
euc.dist <- function(biom, x) {
  xmat <- do.call("rbind", replicate(nrow(biom), x, simplify = FALSE))
  sqrt(rowSums((biom - xmat)^2))
}

#' @keywords internal
#' @param biom Data.frame of biomarker values.
#' @param x Numeric vector, a single grid point.
#' @noRd
standardise_biom_names <- function(biom) {
  biom <- as.data.frame(biom)
  if (is.null(names(biom)) || any(names(biom) == "") ||
      all(grepl("^V[0-9]+$", names(biom)))) {
    names(biom) <- paste0("biom", seq_len(ncol(biom)))
  }
  biom
}
