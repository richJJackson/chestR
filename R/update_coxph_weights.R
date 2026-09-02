#' Update a coxph fit with new case weights
#' @keywords internal
#' @noRd
update_coxph_weights <- function(base, weights) {
  env <- environment(base$formula)
  if (is.null(env) || identical(env, emptyenv())) {
    env <- parent.frame()
  }

  cl <- base$call
  cl$weights <- weights

  if (!is.null(cl$data)) {
    cl$data <- eval(cl$data, envir = env)
  }

  eval(cl, envir = env)
}
