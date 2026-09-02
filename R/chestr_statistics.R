#' Global test statistics from a chestr fit
#'
#' @param x A `"chestr"` object from [chestr()].
#' @param treat_term Treatment coefficient name, or `NULL` to resolve.
#' @param reliable_only If `TRUE`, only reliable grid points enter.
#' @return A list with `T_L2`, `T_MAX`, `beta_global`, `n_points`, `treat_term`.
#' @keywords internal
#' @noRd
chestr_statistics <- function(x, treat_term = NULL, reliable_only = TRUE) {
  if (!inherits(x, "chestr")) {
    stop("x must be a 'chestr' object from chestr().", call. = FALSE)
  }
  treat_term <- resolve_treat_term(x, treat_term)

  se_col <- paste0(treat_term, ".se")
  if (!se_col %in% names(x$estimates)) {
    stop("Standard error column '", se_col, "' not found in x$estimates.",
         call. = FALSE)
  }

  est <- x$estimates
  if (isTRUE(reliable_only) && "reliable" %in% names(est)) {
    est <- est[!is.na(est$reliable) & est$reliable, , drop = FALSE]
  }
  keep <- is.finite(est[[treat_term]]) & is.finite(est[[se_col]]) &
    est[[se_col]] > 0
  est <- est[keep, , drop = FALSE]
  if (nrow(est) == 0L) {
    stop("No finite local estimates available for statistics.", call. = FALSE)
  }

  beta_global <- unname(stats::coef(x$base)[treat_term])
  beta_local <- est[[treat_term]]
  se_local <- pmax(est[[se_col]], 1e-8)

  w <- est$ess_events
  if (is.null(w) || !any(is.finite(w) & w > 0)) {
    w <- est$eff.e
  }
  w[!is.finite(w) | w < 0] <- 0
  if (sum(w) <= 0) w <- rep(1, length(beta_local))
  w <- w / sum(w)

  z_dev <- (beta_local - beta_global) / se_local
  list(
    treat_term = treat_term,
    beta_global = beta_global,
    T_L2 = sum(w * (beta_local - beta_global)^2, na.rm = TRUE),
    T_MAX = max(abs(z_dev), na.rm = TRUE),
    n_points = nrow(est)
  )
}
