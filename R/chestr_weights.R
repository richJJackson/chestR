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
