#' @export
#' @noRd
as.data.frame.chestr <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x$estimates, row.names = row.names, optional = optional, ...)
}
