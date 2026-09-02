#' @export
#' @noRd
print.chestr <- function(x, ...) {
  cat("chestr object\n")
  cat("  biomarkers :", paste(names(x$biom), collapse = ", "), "\n")
  cat("  method     :", x$method, "\n")
  if (!is.null(x$treat_term)) {
    cat("  treat_term :", x$treat_term, "\n")
  }
  cat("  data stored:", if (!is.null(x$data)) "yes" else "no", "\n")
  cat("  grid points:", nrow(x$estimates),
      " (reliable ", sum(x$estimates$reliable, na.rm = TRUE), ")\n", sep = "")
  cat("\nEstimates (head):\n")
  print(utils::head(x$estimates), ...)
  invisible(x)
}
