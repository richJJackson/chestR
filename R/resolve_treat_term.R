#' Resolve treatment coefficient name from a chestr object
#' @keywords internal
#' @noRd
resolve_treat_term <- function(x, treat_term = NULL) {
  if (!is.null(treat_term) && !identical(treat_term, "")) {
    if (!treat_term %in% names(x$estimates)) {
      stop("treat_term '", treat_term, "' not found in x$estimates.",
           call. = FALSE)
    }
    return(treat_term)
  }
  if (!is.null(x$treat_term) && x$treat_term %in% names(x$estimates)) {
    return(x$treat_term)
  }

  coefs <- names(stats::coef(x$base))
  if (length(coefs) == 1L) return(coefs[[1]])

  hit <- grep("treat|trt|arm|vacc|fluvacc|intervention", coefs,
              ignore.case = TRUE, value = TRUE)
  if (length(hit) == 1L) return(hit[[1]])

  stop(
    "Could not infer treat_term. Supply treat_term= or set it in chestr(..., treat_term=).",
    call. = FALSE
  )
}
