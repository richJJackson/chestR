#' Kish effective sample size for a weight vector
#' @keywords internal
#' @noRd
kish_ess <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) return(0)
  sum(w)^2 / sum(w^2)
}
