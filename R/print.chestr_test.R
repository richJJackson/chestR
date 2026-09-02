#' @export
#' @noRd
print.chestr_test <- function(x, ...) {
  cat("chestr permutation test for treatment-effect heterogeneity\n")
  cat("  H0: treatment effect constant over biomarker space\n")
  cat("  treat_term :", x$treat_term, "\n")
  cat("  treat_var  :", x$treat_var, "\n")
  cat("  grid points used:", x$n_points, "\n")
  cat("  permutations    :", x$n_perm_used, "/", x$n_perm_requested, "\n\n")
  tab <- data.frame(
    statistic = c("T_L2", "T_MAX"),
    observed = c(x$T_L2, x$T_MAX),
    p_value = c(x$p_L2, x$p_MAX)
  )
  print(tab, row.names = FALSE, digits = 4, ...)
  invisible(x)
}
