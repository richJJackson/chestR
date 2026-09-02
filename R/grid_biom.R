#' Build a biomarker evaluation grid
#' @keywords internal
#' @noRd
grid_biom <- function(biom, grid.size = 25) {
  dbiom <- as.data.frame(biom)
  grid.cut <- lapply(seq_len(ncol(dbiom)), function(j) {
    seq(min(dbiom[[j]], na.rm = TRUE), max(dbiom[[j]], na.rm = TRUE),
        length.out = grid.size)
  })
  names(grid.cut) <- names(dbiom)
  as.data.frame(expand.grid(grid.cut))
}
