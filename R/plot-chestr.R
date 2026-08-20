#' Plot local treatment effects from a `chestr` object
#'
#' Visualises a local coefficient stored in [chestr()] output using
#' **ggplot2**. Point size reflects local information (effective events);
#' colour reflects the selected coefficient. Model and biomarker data are
#' taken from `x` itself.
#'
#' @param x An object of class `"chestr"` from [chestr()].
#' @param trt.param Column name in `x$estimates` for the effect to colour
#'   (e.g. `"trt"` or `"treat_fGEM"`).
#' @param pnt.scale Multiplier for point size (proportional to local information).
#'   If `NULL`, size is scaled automatically from grid size.
#' @param col.scale Colour limits: `NULL` (symmetric -3 to 3), `"obs"`
#'   (data-driven), or a length-2 numeric range `(min, max)`.
#' @param pts If `TRUE`, overlay observed biomarker values.
#' @param data_cex Point size for observed biomarker data (when `pts = TRUE`).
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
#' \dontrun{
#' plot(cr, trt.param = "trt")
#' }
#'
#' @seealso [chestr()]
plot.chestr <- function(x, trt.param,
                        pnt.scale = 3, col.scale = NULL, pts = TRUE,
                        data_cex = 0.6,
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

  biom <- as.data.frame(x$biom)
  base <- x$base
  n_biom <- ncol(biom)
  if (n_biom < 1L || n_biom > 2L) {
    stop("plot.chestr() supports 1 or 2 biomarker dimensions.", call. = FALSE)
  }

  if (missing(trt.param) || is.null(trt.param)) {
    stop("trt.param must be supplied (column name in x$estimates).", call. = FALSE)
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
    col.lim.lo <- -3
    col.lim.hi <- 3
  } else {
    stop("col.scale must be NULL, 'obs', or a length-2 numeric vector.", call. = FALSE)
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
          size = data_cex, colour = "grey40", alpha = 0.8
        )
    }
    p <- p +
      ggplot2::scale_colour_distiller(
        palette = "RdBu",
        limits = c(col.lim.lo, col.lim.hi),
        oob = scales::squish,
        name = trt.param
      ) +
      ggplot2::scale_size_identity() +
      ggplot2::theme_bw()
  } else {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(.data[[xvar]], .data[[trt.param]])
    ) +
      ggplot2::geom_point(
        ggplot2::aes(
          colour = .data[[trt.param]],
          size = 1.5 * .data[["ps"]]
        ),
        shape = 15, alpha = 0.6
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
      ggplot2::scale_colour_distiller(
        palette = "RdBu",
        limits = c(col.lim.lo, col.lim.hi),
        oob = scales::squish,
        name = trt.param
      ) +
      ggplot2::scale_size_identity() +
      ggplot2::labs(x = xvar, y = trt.param) +
      ggplot2::theme_bw()
  }

  print(p)
  invisible(p)
}

#' @rdname plot.chestr
#' @export
plot_chestr <- function(x, ...) {
  plot.chestr(x, ...)
}
