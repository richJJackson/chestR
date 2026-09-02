#' Permutation test for treatment-effect heterogeneity
#'
#' Tests the null that the treatment effect is constant over biomarker space
#' by permuting treatment labels, re-fitting the global Cox model and
#' [chestr()] surface, and comparing observed heterogeneity statistics to the
#' permutation null.
#'
#' Where possible, arguments are taken from `x`:
#' * `treat_term` from `x$treat_term` (or inferred)
#' * `data` from `x$data` (stored by [chestr()] when recoverable)
#' * `treat_var` inferred from `treat_term`
#' * grid / kernel / ESS settings from `x`
#'
#' Observed statistics are the ESS-weighted L2 deviation and the maximum
#' standardised deviation of local treatment coefficients from the global
#' coefficient.
#'
#' @param x A `"chestr"` object from [chestr()].
#' @param treat_term Treatment coefficient name. Optional if stored on `x` or
#'   uniquely inferable.
#' @param data Analysis data frame. Optional if `x$data` is present or the
#'   data can be recovered from `x$base$call`.
#' @param treat_var Treatment column in `data`. Optional if inferable.
#' @param B Number of permutations (default 100).
#' @param seed Optional random seed for reproducibility.
#' @param reliable_only If `TRUE` (default), only reliable grid points enter
#'   the statistics, and the same ESS rule as `x` is used when rebuilding each
#'   permuted [chestr()] fit.
#'
#' @return An object of class `"chestr_test"` with observed statistics,
#'   permutation p-values, and null draws.
#'
#' @examples
#' \donttest{
#' library(survival)
#' set.seed(1)
#' n <- 80
#' dat <- data.frame(
#'   time = rexp(n, 0.15), status = 1L,
#'   trt = rbinom(n, 1, 0.5),
#'   b1 = rnorm(n), b2 = rnorm(n)
#' )
#' dat$status <- as.integer(dat$time < 4)
#' dat$time[dat$status == 0L] <- 4
#' base <- coxph(Surv(time, status) ~ trt, data = dat)
#' cr <- chestr(base, dat[, c("b1", "b2")], grid.size = 5,
#'              min_events_per_df = 0, treat_term = "trt")
#' tst <- chestr_test(cr, B = 19, seed = 1)
#' tst
#' }
#'
#' @export
#' @seealso [chestr()], [plot.chestr()]
chestr_test <- function(x, treat_term = NULL, data = NULL, treat_var = NULL,
                        B = 100, seed = NULL, reliable_only = TRUE) {
  if (!inherits(x, "chestr")) {
    stop("x must be a 'chestr' object from chestr().", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1) {
    stop("B must be a positive integer.", call. = FALSE)
  }
  B <- as.integer(B)

  treat_term <- resolve_treat_term(x, treat_term)
  obs <- chestr_statistics(x, treat_term = treat_term,
                           reliable_only = reliable_only)

  data <- resolve_chestr_data(x, data)
  treat_var <- resolve_treat_var(x$base, treat_term, treat_var, data)
  if (!treat_var %in% names(data)) {
    stop("treat_var '", treat_var, "' not found in data.", call. = FALSE)
  }

  if (nrow(data) != nrow(x$biom)) {
    stop("nrow(data) must equal nrow(x$biom) (same analysis rows).",
         call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  perm_l2 <- rep(NA_real_, B)
  perm_max <- rep(NA_real_, B)
  fail_msg <- character(0)

  for (b in seq_len(B)) {
    perm_dat <- data
    perm_dat[[treat_var]] <- sample(perm_dat[[treat_var]],
                                    size = nrow(perm_dat),
                                    replace = FALSE)

    res <- tryCatch({
      # Embed the data frame in the call (not the symbol `perm_dat`) so later
      # weight updates can re-eval the fit outside this local frame.
      base_b <- stats::update(x$base, data = perm_dat)
      base_b$call$data <- perm_dat
      cr_b <- chestr(
        base = base_b,
        biom = x$biom,
        grid.size = x$grid.size,
        method = x$method,
        kern.adj = x$kern.adj,
        min_events_per_df = if (isTRUE(reliable_only)) x$min_events_per_df else 0,
        treat_term = treat_term
      )
      chestr_statistics(cr_b, treat_term = treat_term,
                        reliable_only = reliable_only)
    }, error = function(e) e)

    if (inherits(res, "error")) {
      fail_msg <- c(fail_msg, conditionMessage(res))
      next
    }
    perm_l2[b] <- res$T_L2
    perm_max[b] <- res$T_MAX
  }

  valid <- is.finite(perm_l2) & is.finite(perm_max)
  if (!any(valid)) {
    detail <- if (length(fail_msg)) {
      paste0(" Last error: ", fail_msg[length(fail_msg)])
    } else {
      ""
    }
    stop(
      "All permutations failed. Check model stability / ESS settings.",
      detail,
      call. = FALSE
    )
  }
  perm_l2 <- perm_l2[valid]
  perm_max <- perm_max[valid]
  n_used <- length(perm_l2)

  p_l2 <- (1 + sum(perm_l2 >= obs$T_L2)) / (1 + n_used)
  p_max <- (1 + sum(perm_max >= obs$T_MAX)) / (1 + n_used)

  structure(
    list(
      treat_term = treat_term,
      treat_var = treat_var,
      beta_global = obs$beta_global,
      n_points = obs$n_points,
      T_L2 = obs$T_L2,
      T_MAX = obs$T_MAX,
      p_L2 = p_l2,
      p_MAX = p_max,
      null_T_L2 = perm_l2,
      null_T_MAX = perm_max,
      n_perm_requested = B,
      n_perm_used = n_used,
      reliable_only = reliable_only,
      call = match.call()
    ),
    class = "chestr_test"
  )
}
