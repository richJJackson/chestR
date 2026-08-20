#' Global test statistics from a chestr fit
#'
#' Computes the weighted L2 and maximum standardised deviation of local
#' treatment coefficients from the global treatment coefficient.
#'
#' @param x A `"chestr"` object from [chestr()].
#' @param treat_term Column name in `x$estimates` for the treatment effect
#'   (e.g. `"treat_fGEM"` or `"trt"`). If `NULL`, uses `x$treat_term` or an
#'   inferred single treatment coefficient.
#' @param reliable_only If `TRUE` (default), only reliable grid points enter
#'   the statistics.
#'
#' @return A list with `T_L2`, `T_MAX`, `beta_global`, `n_points`, and
#'   `treat_term`.
#'
#' @details
#' Let \eqn{\hat\beta_A} be the global treatment coefficient and
#' \eqn{\hat\beta_A(g)} the local coefficient at grid point \eqn{g}. Define
#' \deqn{T_{L2} = \sum_g \omega(g)\,\{\hat\beta_A(g)-\hat\beta_A\}^2}{T_L2 = sum_g w(g) * (beta(g) - beta)^2}
#' with weights \eqn{\omega(g)} proportional to Kish ESS events at \eqn{g}, and
#' \deqn{T_{MAX} = \max_g |z(g)|}{T_MAX = max |z(g)|}
#' where \eqn{z(g) = \{\hat\beta_A(g)-\hat\beta_A\} / \mathrm{SE}\{\hat\beta_A(g)\}}.
#'
#' @export
#' @seealso [chestr_test()], [chestr()]
chestr_statistics <- function(x, treat_term = NULL, reliable_only = TRUE) {
  if (!inherits(x, "chestr")) {
    stop("x must be a 'chestr' object from chestr().", call. = FALSE)
  }
  treat_term <- resolve_treat_term(x, treat_term)

  se_col <- paste0(treat_term, ".se")
  if (!se_col %in% names(x$estimates)) {
    stop("Standard error column '", se_col, "' not found in x$estimates.",
         call. = FALSE)
  }

  est <- x$estimates
  if (isTRUE(reliable_only) && "reliable" %in% names(est)) {
    est <- est[!is.na(est$reliable) & est$reliable, , drop = FALSE]
  }
  keep <- is.finite(est[[treat_term]]) & is.finite(est[[se_col]]) &
    est[[se_col]] > 0
  est <- est[keep, , drop = FALSE]
  if (nrow(est) == 0L) {
    stop("No finite local estimates available for statistics.", call. = FALSE)
  }

  beta_global <- unname(stats::coef(x$base)[treat_term])
  beta_local <- est[[treat_term]]
  se_local <- pmax(est[[se_col]], 1e-8)

  w <- est$ess_events
  if (is.null(w) || !any(is.finite(w) & w > 0)) {
    w <- est$eff.e
  }
  w[!is.finite(w) | w < 0] <- 0
  if (sum(w) <= 0) w <- rep(1, length(beta_local))
  w <- w / sum(w)

  z_dev <- (beta_local - beta_global) / se_local
  list(
    treat_term = treat_term,
    beta_global = beta_global,
    T_L2 = sum(w * (beta_local - beta_global)^2, na.rm = TRUE),
    T_MAX = max(abs(z_dev), na.rm = TRUE),
    n_points = nrow(est)
  )
}


#' Permutation test for treatment-effect heterogeneity
#'
#' Tests the null that the treatment effect is constant over biomarker space
#' by permuting treatment labels, re-fitting the global Cox model and
#' [chestr()] surface, and comparing observed [chestr_statistics()] values to
#' the permutation null.
#'
#' Where possible, arguments are taken from `x`:
#' * `treat_term` from `x$treat_term` (or inferred)
#' * `data` from `x$data` (stored by [chestr()] when recoverable)
#' * `treat_var` inferred from `treat_term`
#' * grid / kernel / ESS settings from `x`
#'
#' @param x A `"chestr"` object from [chestr()].
#' @param treat_term Treatment coefficient name. Optional if stored on `x` or
#'   uniquely inferable.
#' @param data Analysis data frame. Optional if `x$data` is present or the
#'   data can be recovered from `x$base$call`.
#' @param treat_var Treatment column in `data`. Optional if inferable.
#' @param B Number of permutations (default 199).
#' @param seed Optional random seed for reproducibility.
#' @param reliable_only Passed through to [chestr_statistics()] and used when
#'   rebuilding each permuted [chestr()] fit (same ESS rule as `x`).
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
#' # Uses cr$data, cr$treat_term, and fit settings from cr:
#' tst <- chestr_test(cr, B = 19, seed = 1)
#' tst
#' }
#'
#' @export
#' @seealso [chestr_statistics()], [chestr()]
chestr_test <- function(x, treat_term = NULL, data = NULL, treat_var = NULL,
                        B = 199, seed = NULL, reliable_only = TRUE) {
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


#' @export
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


#' Resolve treatment coefficient name from a chestr object
#' @keywords internal
#' @noRd
resolve_treat_term <- function(x, treat_term = NULL) {
  if (!is.null(treat_term) && !identical(treat_term, "")) {
    if (!treat_term %in% names(x$estimates)) {
      stop("treat_term '", treat_term, "' not found in x$estimates.",
           call. = FALSE)
    }
    return(treat_term)
  }
  if (!is.null(x$treat_term) && x$treat_term %in% names(x$estimates)) {
    return(x$treat_term)
  }

  coefs <- names(stats::coef(x$base))
  if (length(coefs) == 1L) return(coefs[[1]])

  hit <- grep("treat|trt|arm|vacc|fluvacc|intervention", coefs,
              ignore.case = TRUE, value = TRUE)
  if (length(hit) == 1L) return(hit[[1]])

  stop(
    "Could not infer treat_term. Supply treat_term= or set it in chestr(..., treat_term=).",
    call. = FALSE
  )
}


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


#' Infer treatment variable name from a coefficient name
#' @keywords internal
#' @noRd
resolve_treat_var <- function(base, treat_term, treat_var, data) {
  if (!is.null(treat_var)) return(treat_var)

  if (treat_term %in% names(data)) return(treat_term)

  mf <- tryCatch(stats::model.frame(base), error = function(e) NULL)
  nms <- if (!is.null(mf)) names(mf) else names(data)
  if (!is.null(mf) && length(nms) > 0L) {
    nms <- nms[-1L]
  }

  for (nm in nms) {
    v <- if (!is.null(mf) && nm %in% names(mf)) mf[[nm]] else data[[nm]]
    if (is.null(v)) next
    if (is.factor(v) || is.character(v)) {
      lev <- if (is.factor(v)) levels(v) else unique(as.character(v))
      if (treat_term %in% paste0(nm, lev)) return(nm)
    }
  }

  hits <- nms[startsWith(treat_term, nms)]
  if (length(hits) >= 1L) {
    return(hits[which.max(nchar(hits))])
  }

  stop(
    "Could not infer treat_var from treat_term '", treat_term,
    "'. Please supply treat_var=.",
    call. = FALSE
  )
}
