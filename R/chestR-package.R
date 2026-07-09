#' chestR: Kernel-Weighted Cox Regression
#'
#' @description
#' The **chestR** package implements kernel-weighted Cox regression for
#' exploring treatment effect heterogeneity and identifying candidate predictive
#' biomarkers. A global Cox model is re-fitted at each point on a biomarker grid
#' using Gaussian (or inverse-square) kernel weights; local coefficient estimates
#' are then visualised across biomarker space.
#'
#' @details
#' The main workflow is:
#'
#' 1. Fit a global `coxph` model (`base`) including treatment and covariates.
#' 2. Call [chestr()] with biomarker values to obtain local estimates on a grid.
#' 3. Use [plot_chestr()] to visualise how a coefficient (e.g. treatment) varies.
#'
#' @seealso [chestr()], [plot_chestr()], [grid_biom()]
#' @importFrom stats dnorm sd update
#' @importFrom survival Surv
#' @keywords internal
"_PACKAGE"
