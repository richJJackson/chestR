#' Plot local treatment effects across biomarker space
#'
#' Visualises a local coefficient from [chestr()] output over one or two
#' biomarker dimensions. Point size reflects local information (effective
#' events); colour reflects the magnitude of the selected coefficient.
#'
#' @param cr.ob Output from [chestr()].
#' @param trt.param Column name in `cr.ob` for the effect to colour points
#'   (e.g. `"ArmGEM"` or `"as.factor.fluvacc.1"`).
#' @param base The fitted `coxph` object passed to [chestr()].
#' @param biom Biomarker values used in [chestr()].
#' @param pnt.scale Multiplier for point size (proportional to local information).
#'   If `NULL`, size is scaled automatically from grid size.
#' @param col.scale Colour breaks: `NULL` (symmetric -3 to 3), `"obs"`
#'   (data-driven), or a length-2 numeric range `(min, max)`.
#' @param pts If `TRUE`, overlay observed biomarker values.
#' @param data_cex Point size for observed biomarker data (when `pts = TRUE`).
#' @param data_col Point colour for observed biomarker data.
#' @param reliable_only If `TRUE` (default), only plot grid points that passed
#'   the events-per-df safeguard.
#'
#' @details
#' Objects returned by [chestr()] have S3 class `"chestr"`, so
#' `plot(cr, trt.param = ..., base = ..., biom = ...)` dispatches to the
#' `plot` method. Direct calls to `plot.chestr()` or [plot_chestr()] are
#' equivalent.
#'
#' @return Invisibly returns `cr.ob`.
#'
#' @export
#' @name plot_chestr
#'
#' @examples
#' \dontrun{
#' plot_chestr(cr, trt.param = "trt", base = base, biom = biom)
#' }
#'
#' @seealso [chestr()]
plot_chestr <- function(cr.ob, trt.param, base, biom,
                        pnt.scale = 3, col.scale = NULL, pts = TRUE,
                        data_cex = 0.45,
                        data_col = grDevices::rgb(0, 0, 0, 0.25),
                        reliable_only = TRUE) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required for plot_chestr().", call. = FALSE)
  }

  biom <- as.data.frame(biom)
  n_biom <- ncol(biom)
  if (n_biom < 1L || n_biom > 2L) {
    stop("plot_chestr() supports 1 or 2 biomarker dimensions.", call. = FALSE)
  }

  if (!trt.param %in% names(cr.ob)) {
    stop("trt.param '", trt.param, "' not found in chestr output.", call. = FALSE)
  }

  plot_df <- cr.ob
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

  invisible(cr.ob)
}


#' @rdname plot_chestr
#' @param x Output from [chestr()] (S3 class `"chestr"`).
#' @param ... Additional arguments passed to [plot_chestr()]
#'   (`trt.param`, `base`, `biom`, and plotting options).
#' @method plot chestr
#' @export
#' @rawNamespace export(plot.chestr)
plot.chestr <- function(x, ...) {
  plot_chestr(x, ...)
}
