#' Infer treatment variable name from a coefficient name
#' @keywords internal
#' @noRd
resolve_treat_var <- function(base, treat_term, treat_var, data) {
  if (!is.null(treat_var)) return(treat_var)

  if (treat_term %in% names(data)) return(treat_term)

  mf <- tryCatch(stats::model.frame(base), error = function(e) NULL)
  nms <- if (!is.null(mf)) names(mf) else names(data)
  if (!is.null(mf) && length(nms) > 0L) {
    nms <- nms[-1L]
  }

  for (nm in nms) {
    v <- if (!is.null(mf) && nm %in% names(mf)) mf[[nm]] else data[[nm]]
    if (is.null(v)) next
    if (is.factor(v) || is.character(v)) {
      lev <- if (is.factor(v)) levels(v) else unique(as.character(v))
      if (treat_term %in% paste0(nm, lev)) return(nm)
    }
  }

  hits <- nms[startsWith(treat_term, nms)]
  if (length(hits) >= 1L) {
    return(hits[which.max(nchar(hits))])
  }

  stop(
    "Could not infer treat_var from treat_term '", treat_term,
    "'. Please supply treat_var=.",
    call. = FALSE
  )
}
