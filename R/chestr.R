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
#'
#' @return A data.frame with grid coordinates, coefficient estimates, standard
#'   errors (`*.se` columns), `eff.n` (effective sample size), and `eff.e` /
#'   `eff.event` (effective weighted events).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data <- readRDS("path/to/data.rds")
#' base <- coxph(Surv(time, status) ~ trt + cov1, data = data)
#' cr <- chestr(base, data$biomarker, grid.size = 20)
#' }
#'
#' @seealso [chestr_point()], [plot_chestr()], [grid_biom()]
chestr <- function(base, biom, grid.size = 25,
                   method = c("distance", "square_distance", "legacy"),
                   kern.adj = 4, bw = NULL) {
  method <- match.arg(method)
  biom <- standardise_biom_names(biom)

  grid.xy <- grid_biom(biom, grid.size = grid.size)

  wmod_ls <- lapply(seq_len(nrow(grid.xy)), function(i) {
    chestr_point(base, biom, grid.xy[i, , drop = FALSE],
                 bw = bw, method = method, kern.adj = kern.adj)
  })

  eff.n <- vapply(wmod_ls, function(x) sum(x$weights), numeric(1))
  eff.e <- vapply(wmod_ls, function(x) sum(x$y[, 2] * x$weights), numeric(1))

  coef_mat <- do.call(rbind, lapply(wmod_ls, function(x) {
    sc <- summary(x)$coef
    vals <- as.data.frame(t(sc[, 1, drop = FALSE]), stringsAsFactors = FALSE)
    names(vals) <- rownames(sc)
    vals
  }))
  se_mat <- do.call(rbind, lapply(wmod_ls, function(x) {
    sc <- summary(x)$coef
    vals <- as.data.frame(t(sc[, 3, drop = FALSE]), stringsAsFactors = FALSE)
    names(vals) <- rownames(sc)
    vals
  }))
  names(se_mat) <- paste(names(se_mat), "se", sep = ".")

  ret <- cbind(grid.xy, coef_mat, se_mat,
               eff.n = eff.n, eff.e = eff.e, eff.event = eff.e)
  rownames(ret) <- NULL
  ret
}

#' Fit a weighted Cox model at a single grid point
#'
#' Internal workhorse for [chestr()]. Exported for advanced use (e.g. evaluating
#' a single biomarker profile).
#'
#' @inheritParams chestr
#' @param x Numeric vector giving a single grid point (one value per biomarker).
#'
#' @return An updated `coxph` object fitted with kernel weights at `x`.
#'
#' @export
#' @seealso [chestr()]
chestr_point <- function(base, biom, x, bw = NULL,
                         method = c("distance", "square_distance", "legacy"),
                         kern.adj = 4) {
  method <- match.arg(method)
  biom <- as.data.frame(biom)
  x <- as.numeric(x)

  if (method == "legacy") {
    dist <- euc.dist(biom, x)
    kern.sd <- mean(vapply(seq_len(ncol(biom)), function(j) {
      sd(biom[[j]], na.rm = TRUE)
    }, numeric(1)), na.rm = TRUE) / kern.adj
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

  update_coxph_weights(base, weight + 1e-5)
}