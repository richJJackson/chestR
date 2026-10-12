#' Plot local treatment effects from a `chestr` object
#'
#' Visualises a local coefficient stored in [chestr()] output using
#' **ggplot2**. For two biomarkers, a coloured grid of local log-hazard
#' ratios is shown (point size reflects local information), with observed
#' patients overlaid as small black points with a white outline when
#' `pts = TRUE`. The default colour scale is blue–white–red and centred on
#' the unweighted treatment coefficient from `x$base`, so white is that
#' coefficient and colour is the deviation around it. For a single
#' biomarker, a line of local hazard
#' ratios is shown with a dashed horizontal reference at the global hazard
#' ratio from `x$base`. Model and biomarker data are taken from `x` itself.
#'
#' @param x An object of class `"chestr"` from [chestr()].
#' @param trt.param Column name in `x$estimates` for the effect to colour
#'   (2D) or plot as a hazard-ratio curve (1D). If `NULL`, uses
#'   `x$treat_term` when set / uniquely inferable.
#' @param pnt.scale Multiplier for point size (proportional to local information).
#'   If `NULL`, size is scaled automatically from grid size. Used for 2D plots.
#' @param col.scale Colour limits for 2D plots. `NULL` (default) uses a
#'   symmetric window around the unweighted `trt.param` coefficient from
#'   `x$base`, wide enough to cover -1 and 1, so white is the legend
#'   midpoint. `"obs"` uses the range of the plotted local estimates. A
#'   length-2 numeric vector sets `(min, max)` directly. Whenever that
#'   coefficient lies inside the limits, it is drawn as white.
#' @param pts If `TRUE`, overlay observed biomarker values (2D points or 1D rug).
#' @param data_cex Point size for observed biomarker data (when `pts = TRUE`).
#'   On a 2D surface these are black points with a white outline.
#' @param reliable_only If `TRUE` (default), only plot grid points that passed
#'   the events-per-df safeguard.
#' @param ... Ignored (for S3 compatibility).
#'
#' @return The ggplot object (invisibly).
#'
#' @export
#' @method plot chestr
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
#' cr <- chestr(base, dat$biom, grid.size = 5, min_events_per_df = 0,
#'              treat_term = "trt")
#' head(cr$estimates)
#' plot(cr, trt.param = "trt")
#' }
#'
#' @seealso [chestr()], [chestr_test()]
plot.chestr <- function(x, trt.param = NULL,
                        pnt.scale = 3, col.scale = NULL, pts = TRUE,
                        data_cex = 1.6,
                        reliable_only = TRUE,
                        ...) {
  if (!inherits(x, "chestr")) {
    stop("x must be a 'chestr' object from chestr().", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.chestr().", call. = FALSE)
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required for plot.chestr().", call. = FALSE)
  }

  if (is.null(trt.param)) {
    trt.param <- tryCatch(resolve_treat_term(x, NULL), error = function(e) NULL)
  }
  if (is.null(trt.param)) {
    stop("trt.param must be supplied (or set treat_term in chestr()).",
         call. = FALSE)
  }

  biom <- as.data.frame(x$biom)
  base <- x$base
  n_biom <- ncol(biom)
  if (n_biom < 1L || n_biom > 2L) {
    stop("plot.chestr() supports 1 or 2 biomarker dimensions.", call. = FALSE)
  }

  if (!trt.param %in% names(x$estimates)) {
    stop("trt.param '", trt.param, "' not found in x$estimates.", call. = FALSE)
  }

  plot_df <- x$estimates
  if (isTRUE(reliable_only) && "reliable" %in% names(plot_df)) {
    plot_df <- plot_df[!is.na(plot_df$reliable) & plot_df$reliable, , drop = FALSE]
  }
  plot_df <- plot_df[is.finite(plot_df[[trt.param]]), , drop = FALSE]
  if (nrow(plot_df) == 0L) {
    stop("No reliable/finite local estimates available to plot.", call. = FALSE)
  }

  nev <- sum(base$y[, 2])
  plot_df$ps <- plot_df$eff.e / nev
  if (is.null(pnt.scale)) {
    scl <- sqrt(nrow(plot_df)) / 15
    plot_df$ps <- plot_df$ps * scl / max(plot_df$ps)
  } else {
    plot_df$ps <- plot_df$ps * pnt.scale / max(plot_df$ps)
  }

  out <- plot_df[[trt.param]]
  beta_global <- suppressWarnings(unname(stats::coef(base)[trt.param]))
  if (identical(col.scale, "obs")) {
    col.lim.lo <- min(out)
    col.lim.hi <- max(out)
  } else if (is.numeric(col.scale)) {
    if (length(col.scale) != 2L) {
      stop("col.scale must be length 2 (min and max) when numeric.", call. = FALSE)
    }
    col.lim.lo <- col.scale[1]
    col.lim.hi <- col.scale[2]
  } else if (is.null(col.scale)) {
    # Symmetric about the unweighted coefficient, and wide enough to
    # include -1 and 1, so white sits at that coefficient.
    if (!is.finite(beta_global)) {
      if (n_biom == 1L) {
        col.lim.lo <- -1
        col.lim.hi <- 1
      } else {
        stop("Global coefficient for '", trt.param, "' is not available.",
             call. = FALSE)
      }
    } else {
      half <- max(abs(beta_global + 1), abs(beta_global - 1))
      col.lim.lo <- beta_global - half
      col.lim.hi <- beta_global + half
    }
  } else {
    stop("col.scale must be NULL, 'obs', or a length-2 numeric vector.", call. = FALSE)
  }
  if (is.finite(beta_global) &&
      beta_global >= col.lim.lo && beta_global <= col.lim.hi) {
    col_mid <- beta_global
  } else {
    col_mid <- (col.lim.lo + col.lim.hi) / 2
  }

  xvar <- names(biom)[1]

  if (n_biom == 2L) {
    yvar <- names(biom)[2]
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = plot_df,
        ggplot2::aes(
          .data[[xvar]], .data[[yvar]],
          colour = .data[[trt.param]], size = 1.5 * .data[["ps"]]
        ),
        shape = 15, alpha = 0.6
      )
    if (isTRUE(pts)) {
      p <- p +
        ggplot2::geom_point(
          data = biom,
          ggplot2::aes(.data[[xvar]], .data[[yvar]]),
          shape = 21, size = data_cex, stroke = 0.35,
          fill = "black", colour = "white"
        )
    }
    p <- p +
      ggplot2::scale_colour_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = col_mid,
        limits = c(col.lim.lo, col.lim.hi),
        breaks = if (is.null(col.scale)) {
          c(col.lim.lo, col_mid, col.lim.hi)
        } else {
          ggplot2::waiver()
        },
        labels = if (is.null(col.scale)) {
          function(br) sprintf("%.3f", br)
        } else {
          ggplot2::waiver()
        },
        oob = scales::squish,
        name = trt.param
      ) +
      ggplot2::scale_size_identity() +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "grey92")
      )
  } else {
    plot_df <- plot_df[order(plot_df[[xvar]]), , drop = FALSE]
    if (!is.finite(beta_global)) {
      stop("Global coefficient for '", trt.param, "' is not available.",
           call. = FALSE)
    }
    hr_global <- exp(beta_global)
    plot_df$.local_hr <- exp(plot_df[[trt.param]])

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(.data[[xvar]], .data[[".local_hr"]])
    ) +
      ggplot2::geom_line(linewidth = 1, colour = "#2166AC") +
      ggplot2::geom_hline(
        yintercept = hr_global,
        linetype = "dashed",
        colour = "grey30",
        linewidth = 0.6
      ) +
      ggplot2::annotate(
        "text",
        x = max(plot_df[[xvar]], na.rm = TRUE),
        y = hr_global,
        label = sprintf("Global HR = %.2f", hr_global),
        hjust = 1,
        vjust = -0.6,
        size = 3.2,
        colour = "grey30"
      )
    if (isTRUE(pts)) {
      p <- p +
        ggplot2::geom_rug(
          data = biom,
          ggplot2::aes(x = .data[[xvar]]),
          inherit.aes = FALSE,
          sides = "b",
          alpha = 0.5,
          colour = "grey40"
        )
    }
    p <- p +
      ggplot2::labs(x = xvar, y = "Hazard ratio") +
      ggplot2::theme_bw()
  }

  print(p)
  invisible(p)
}
