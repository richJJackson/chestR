#' Recover analysis data from a coxph fit when row count matches biom
#' @keywords internal
#' @noRd
extract_base_data <- function(base, n_biom) {
  cl <- base$call
  if (is.null(cl$data)) return(NULL)
  env <- environment(stats::formula(base))
  if (is.null(env) || identical(env, emptyenv())) {
    env <- parent.frame(2)
  }
  dat <- tryCatch(eval(cl$data, envir = env), error = function(e) NULL)
  if (is.null(dat)) return(NULL)
  dat <- as.data.frame(dat)
  if (nrow(dat) != n_biom) return(NULL)
  dat
}
