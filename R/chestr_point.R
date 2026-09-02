#' Fit a weighted Cox model at a single grid point
#'
#' @inheritParams chestr
#' @param x Numeric vector giving a single grid point (one value per biomarker).
#' @param warn If `TRUE`, warn when the events-per-df safeguard skips the fit.
#' @return A list with `fit`, `eff.n`, `eff.e`, `ess_events`, `events_per_df`,
#'   and `skipped`.
#' @keywords internal
#' @noRd
chestr_point <- function(base, biom, x, bw = NULL,
                         method = c("distance", "square_distance", "legacy"),
                         kern.adj = 4,
                         min_events_per_df = 10,
                         warn = TRUE) {
  method <- match.arg(method)
  biom <- as.data.frame(biom)
  x <- as.numeric(x)

  weight <- chestr_weights(biom, x, method = method, kern.adj = kern.adj, bw = bw)
  events <- as.numeric(base$y[, 2])
  eff.n <- sum(weight, na.rm = TRUE)
  eff.e <- sum(events * weight, na.rm = TRUE)
  ess_events <- kish_ess(weight[events == 1])
  n_df <- length(stats::coef(base))
  events_per_df <- ess_events / n_df

  thr <- if (is.null(min_events_per_df)) 0 else min_events_per_df
  if (events_per_df < thr) {
    if (isTRUE(warn)) {
      warning(
        "Skipping local fit: Kish effective events/df = ",
        signif(events_per_df, 3), " < ", thr,
        " (Kish ESS events = ", signif(ess_events, 3),
        ", df = ", n_df, ").",
        call. = FALSE
      )
    }
    return(list(
      fit = NULL,
      eff.n = eff.n,
      eff.e = eff.e,
      ess_events = ess_events,
      events_per_df = events_per_df,
      skipped = TRUE
    ))
  }

  fit <- update_coxph_weights(base, weight + 1e-5)
  list(
    fit = fit,
    eff.n = eff.n,
    eff.e = eff.e,
    ess_events = ess_events,
    events_per_df = events_per_df,
    skipped = FALSE
  )
}
