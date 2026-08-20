#' Plot local treatment effects from a `chestr` object
#'
#' Visualises a local coefficient stored in [chestr()] output. Point size
#' reflects local information (effective events); colour reflects the selected
#' coefficient. Model and biomarker data are taken from `x` itself.
#'
#' @param x An object of class `"chestr"` from [chestr()].
#' @param trt.param Column name in `x$estimates` for the effect to colour
#'   (e.g. `"trt"` or `"treat_fGEM"`).
#' @param pnt.scale Multiplier for point size (proportional to local information).
#'   If `NULL`, size is scaled automatically from grid size.
#' @param col.scale Colour breaks: `NULL` (symmetric -3 to 3), `"obs"`
#'   (data-driven), or a length-2 numeric range `(min, max)`.
#' @param pts If `TRUE`, overlay observed biomarker values.
#' @param data_cex Point size for observed biomarker data (when `pts = TRUE`).
#' @param data_col Point colour for observed biomarker data.
#' @param reliable_only If `TRUE` (default), only plot grid points that passed
#'   the events-per-df safeguard.
#' @param ... Ignored (for S3 compatibility).
#'
#' @return Invisibly returns `x`.
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
                        data_cex = 0.45,
                        data_col = grDevices::rgb(0, 0, 0, 0.25),
                        reliable_only = TRUE,
                        ...) {
  if (!inherits(x, "chestr")) {
    stop("x must be a 'chestr' object from chestr().", call. = FALSE)
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required for plot.chestr().", call. = FALSE)
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
    col.breaks <- c(-Inf, seq(min(out), max(out), length.out = 10), Inf)
  } else if (is.numeric(col.scale)) {
    if (length(col.scale) != 2L) {
      stop("col.scale must be length 2 (min and max) when numeric.", call. = FALSE)
    }
    col.breaks <- c(-Inf, seq(col.scale[1], col.scale[2], length.out = 10), Inf)
  } else if (is.null(col.scale)) {
    col.breaks <- c(-Inf, seq(-3, 3, length.out = 10), Inf)
  } else {
    stop("col.scale must be NULL, 'obs', or a length-2 numeric vector.", call. = FALSE)
  }

  n_bins <- length(col.breaks) - 1L
  pal <- RColorBrewer::brewer.pal(max(3L, n_bins), "RdBu")
  if (length(pal) > n_bins) pal <- pal[seq_len(n_bins)]
  cl.lev <- cut(out, breaks = col.breaks)
  cl.cls <- cut(out, breaks = col.breaks, labels = pal)

  if (n_biom == 2L) {
    grid <- plot_df[, seq_len(2), drop = FALSE]
    graphics::par(mar = c(4, 4, 5, 15))
    if (pts) {
      graphics::plot(biom, type = "p", pch = 16, cex = data_cex, col = data_col)
    } else {
      graphics::plot(biom, type = "n")
    }
    graphics::points(grid, cex = plot_df$ps, pch = 15, col = as.character(cl.cls))
    mc <- apply(grid, 2, max, na.rm = TRUE)
    graphics::legend(mc[1] * 1.1, mc[2], levels(cl.lev), col = levels(cl.cls),
                     pch = 20, xpd = TRUE, bty = "n")
  } else {
    grid_x <- plot_df[[1]]
    graphics::par(mar = c(4, 4, 5, 12))
    graphics::plot(grid_x, out, type = "n",
                   xlab = names(biom)[1], ylab = trt.param)
    graphics::points(grid_x, out, cex = plot_df$ps, pch = 15,
                     col = as.character(cl.cls))
    if (pts) graphics::rug(biom[[1]])
    graphics::legend(x = max(grid_x) * 1.05, y = max(out, na.rm = TRUE),
                     legend = levels(cl.lev), col = levels(cl.cls),
                     pch = 20, xpd = TRUE, bty = "n")
  }

  invisible(x)
}

#' @rdname plot.chestr
#' @export
plot_chestr <- function(x, ...) {
  plot.chestr(x, ...)
}
