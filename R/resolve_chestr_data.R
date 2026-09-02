#' Resolve analysis data for permutation tests
#' @keywords internal
#' @noRd
resolve_chestr_data <- function(x, data) {
  if (!is.null(data)) {
    return(as.data.frame(data))
  }
  if (!is.null(x$data)) {
    return(as.data.frame(x$data))
  }
  cl <- x$base$call
  if (is.null(cl$data)) {
    stop("data not available on x and x$base$call has no data= argument. ",
         "Supply data= to chestr_test().",
         call. = FALSE)
  }
  env <- environment(stats::formula(x$base))
  if (is.null(env) || identical(env, emptyenv())) {
    env <- parent.frame(2)
  }
  dat <- tryCatch(eval(cl$data, envir = env), error = function(e) NULL)
  if (is.null(dat)) {
    stop("Could not recover data from x$base$call; please supply data=.",
         call. = FALSE)
  }
  as.data.frame(dat)
}
