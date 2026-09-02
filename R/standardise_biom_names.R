#' Ensure biomarker columns have usable names
#' @keywords internal
#' @noRd
standardise_biom_names <- function(biom) {
  biom <- as.data.frame(biom)
  if (is.null(names(biom)) || any(names(biom) == "") ||
      all(grepl("^V[0-9]+$", names(biom)))) {
    names(biom) <- paste0("biom", seq_len(ncol(biom)))
  }
  biom
}
