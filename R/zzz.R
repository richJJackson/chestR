#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  assign("plot.chestr", plot_chestr, envir = asNamespace(pkgname))
}
