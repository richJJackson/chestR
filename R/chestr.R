#' Kernel-weighted Cox estimates over a biomarker grid
#'
#' Fits a weighted Cox model at each point on a biomarker grid. Observations
#' near the grid point receive higher weight, yielding local estimates of
#' treatment and covariate effects.
#'
#' @param base A fitted `coxph` object (the global model).
#' @param biom Biomarker values: a vector, matrix, or data.frame (1--2 columns
#'   supported for [plot_chestr()]).
#' @param grid.size Number of grid points per biomarker dimension.
#' @param method Weighting scheme:
#'   * `"distance"` (default): scaled biomarkers + Gaussian kernel.
#'   * `"square_distance"`: inverse-square distance weights.
#'   * `"legacy"`: unscaled Gaussian kernel; use `kern.adj` for bandwidth.
#' @param kern.adj Bandwidth divisor when `method = "legacy"`. Larger values
#'   give narrower kernels.
#' @param bw Optional kernel standard deviation; reserved for future use.
#' @param min_events_per_df Minimum Kish effective events per model degree of
#'   freedom required to fit a local model (default `10`). Effective events use
#'   Kish ESS among event weights: `(sum w_e)^2 / sum(w_e^2)`. Grid points below
#'   the threshold are skipped (coefficients set to `NA`) and a warning is
#'   issued. Set to `0` or `NULL` to disable.
#'
#' @return A data.frame with grid coordinates, coefficient estimates, standard
#'   errors (`*.se` columns), `eff.n`, `eff.e`/`eff.event` (sum of event
#'   weights), `ess_events` (Kish ESS), `events_per_df`, and `reliable`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(survival)
#' set.seed(1)
#' n <- 60
#' dat <- data.frame(
#'   time = rexp(n, 0.1), status = 1L,
#'   trt = rbinom(n, 1, 0.5), biom = rnorm(n)
#' )
#' dat$status <- as.integer(dat$time < 5)
#' dat$time[dat$status == 0L] <- 5
#' base <- coxph(Surv(time, status) ~ trt, data = dat)
#' cr <- chestr(base, dat$biom, grid.size = 5, min_events_per_df = 0)
#' head(cr)
#' }
#'
#' @seealso [chestr_point()], [plot_chestr()], [grid_biom()]
chestr <- function(base, biom, grid.size = 25,
                   method = c("distance", "square_distance", "legacy"),
                   kern.adj = 4, bw = NULL,
                   min_events_per_df = 10) {
  method <- match.arg(method)
  biom <- standardise_biom_names(biom)

  if (!is.null(min_events_per_df) &&
      (!is.numeric(min_events_per_df) || length(min_events_per_df) != 1L ||
       is.na(min_events_per_df) || min_events_per_df < 0)) {
    stop("min_events_per_df must be NULL or a single non-negative number.",
         call. = FALSE)
  }

  n_df <- length(stats::coef(base))
  if (n_df < 1L) {
    stop("base model has no estimated coefficients.", call. = FALSE)
  }

  grid.xy <- grid_biom(biom, grid.size = grid.size)
  term_names <- names(stats::coef(base))

  point_res <- lapply(seq_len(nrow(grid.xy)), function(i) {
    chestr_point(base, biom, grid.xy[i, , drop = FALSE],
                 bw = bw, method = method, kern.adj = kern.adj,
                 min_events_per_df = min_events_per_df,
                 warn = FALSE)
  })

  skipped <- vapply(point_res, function(x) isTRUE(x$skipped), logical(1))
  n_skipped <- sum(skipped)
  if (n_skipped > 0L) {
    thr <- if (is.null(min_events_per_df)) 0 else min_events_per_df
    warning(
      "Skipped ", n_skipped, " of ", length(point_res),
      " grid point(s) with fewer than ", thr,
      " Kish effective events per model df (", n_df, " df; require >= ",
      thr * n_df, " Kish ESS events). Coefficients set to NA.",
      call. = FALSE
    )
  }
  if (n_skipped == length(point_res)) {
    stop(
      "All grid points were skipped by the events-per-df safeguard. ",
      "Lower min_events_per_df, increase bandwidth, or reduce model complexity.",
      call. = FALSE
    )
  }

  eff.n <- vapply(point_res, function(x) x$eff.n, numeric(1))
  eff.e <- vapply(point_res, function(x) x$eff.e, numeric(1))
  ess_events <- vapply(point_res, function(x) x$ess_events, numeric(1))
  events_per_df <- vapply(point_res, function(x) x$events_per_df, numeric(1))

  coef_mat <- as.data.frame(matrix(NA_real_, nrow = length(point_res),
                                   ncol = length(term_names)))
  se_mat <- as.data.frame(matrix(NA_real_, nrow = length(point_res),
                                 ncol = length(term_names)))
  names(coef_mat) <- term_names
  names(se_mat) <- paste(term_names, "se", sep = ".")

  for (i in seq_along(point_res)) {
    fit <- point_res[[i]]$fit
    if (is.null(fit)) next
    sc <- summary(fit)$coef
    coef_mat[i, ] <- sc[term_names, 1]
    se_mat[i, ] <- sc[term_names, 3]
  }

  ret <- cbind(grid.xy, coef_mat, se_mat,
               eff.n = eff.n, eff.e = eff.e, eff.event = eff.e,
               ess_events = ess_events,
               events_per_df = events_per_df,
               reliable = !skipped)
  rownames(ret) <- NULL
  class(ret) <- c("chestr", class(ret))
  ret
}

#' Fit a weighted Cox model at a single grid point
#'
#' Workhorse for [chestr()]. Exported for advanced use (e.g. evaluating
#' a single biomarker profile).
#'
#' @inheritParams chestr
#' @param x Numeric vector giving a single grid point (one value per biomarker).
#' @param warn If `TRUE`, warn when the events-per-df safeguard skips the fit.
#'
#' @return A list with:
#'   * `fit`: fitted `coxph` object, or `NULL` if skipped
#'   * `eff.n`, `eff.e`, `ess_events`, `events_per_df`
#'   * `skipped`: logical
#'
#' @export
#' @seealso [chestr()]
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
  eff.n <- sum(weight)
  eff.e <- sum(events * weight)
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

#' Kish effective sample size for a weight vector
#' @keywords internal
#' @noRd
kish_ess <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (length(w) == 0L) return(0)
  sum(w)^2 / sum(w^2)
}

#' Compute kernel weights for a single grid point
#' @keywords internal
#' @noRd
chestr_weights <- function(biom, x,
                           method = c("distance", "square_distance", "legacy"),
                           kern.adj = 4, bw = NULL) {
  method <- match.arg(method)
  biom <- as.data.frame(biom)
  x <- as.numeric(x)

  if (method == "legacy") {
    dist <- euc.dist(biom, x)
    kern.sd <- if (!is.null(bw)) {
      bw
    } else {
      mean(vapply(seq_len(ncol(biom)), function(j) {
        stats::sd(biom[[j]], na.rm = TRUE)
      }, numeric(1)), na.rm = TRUE) / kern.adj
    }
    weight <- stats::dnorm(dist, 0, kern.sd) / stats::dnorm(0, 0, kern.sd)
  } else {
    sbiom <- scale(biom)
    center <- attr(sbiom, "scaled:center")
    scale_fac <- attr(sbiom, "scaled:scale")
    x_scaled <- (x - center) / scale_fac
    sbiom <- as.data.frame(sbiom)
    dist <- euc.dist(sbiom, x_scaled)

    weight <- if (method == "square_distance") {
      1 / pmax(dist^2, .Machine$double.eps)
    } else {
      stats::dnorm(dist, 0, 1)
    }
  }
  weight
}
