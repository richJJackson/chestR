test_that("grid_biom returns expected dimensions", {
  g <- grid_biom(data.frame(x = 1:5, y = 11:15), grid.size = 4)
  expect_equal(nrow(g), 16)
  expect_equal(names(g), c("x", "y"))
})

test_that("chestr returns classed list with estimates, base, biom", {
  skip_if_not_installed("survival")
  set.seed(42)
  n <- 80
  dat <- data.frame(
    st = rexp(n, rate = 0.1),
    trt = rbinom(n, 1, 0.5),
    b1 = rnorm(n),
    b2 = rnorm(n)
  )
  dat$cen <- as.numeric(dat$st < 5)
  dat$st[dat$cen == 0] <- 5
  base <- survival::coxph(survival::Surv(st, cen) ~ trt, data = dat)
  biom <- dat[, c("b1", "b2")]
  cr <- suppressWarnings(
    chestr(base, biom, grid.size = 5,
           method = "legacy", kern.adj = 2, min_events_per_df = 5)
  )
  expect_s3_class(cr, "chestr")
  expect_true(is.list(cr))
  expect_true(all(c("estimates", "base", "biom") %in% names(cr)))
  expect_identical(cr$base, base)
  expect_equal(names(cr$biom), c("b1", "b2"))
  expect_true("trt" %in% names(cr$estimates))
  expect_true("eff.event" %in% names(cr$estimates))
  expect_true("ess_events" %in% names(cr$estimates))
  expect_equal(nrow(cr$estimates), 25)
  expect_true(any(grepl("\\.se$", names(cr$estimates))))
})

test_that("chestr assigns default biomarker names for unnamed inputs", {
  skip_if_not_installed("survival")
  set.seed(1)
  n <- 40
  dat <- data.frame(
    st = rexp(n, 0.1),
    x = rnorm(n),
    b1 = rnorm(n),
    b2 = rnorm(n)
  )
  dat$cen <- 1L
  b <- as.matrix(dat[, c("b1", "b2")])
  dimnames(b) <- NULL
  base <- survival::coxph(survival::Surv(st, cen) ~ x, data = dat)
  cr <- chestr(base, b, grid.size = 3, method = "legacy", min_events_per_df = 0)
  expect_true(all(c("biom1", "biom2") %in% names(cr$estimates)))
  expect_true(all(c("biom1", "biom2") %in% names(cr$biom)))
  expect_true("x" %in% names(cr$estimates))
})

test_that("chestr skips grid points below min_events_per_df", {
  skip_if_not_installed("survival")
  set.seed(99)
  n <- 200
  dat <- data.frame(
    st = rexp(n, rate = 0.15),
    trt = rbinom(n, 1, 0.5),
    b1 = rnorm(n),
    b2 = rnorm(n)
  )
  dat$cen <- as.numeric(dat$st < 5)
  dat$st[dat$cen == 0] <- 5
  base <- survival::coxph(survival::Surv(st, cen) ~ trt, data = dat)

  expect_warning(
    cr <- chestr(base, dat[, c("b1", "b2")], grid.size = 8,
                 method = "legacy", kern.adj = 2, min_events_per_df = 20),
    "Skipped"
  )
  est <- cr$estimates
  expect_true("reliable" %in% names(est))
  expect_true("ess_events" %in% names(est))
  expect_true(any(est$reliable))
  expect_true(any(!est$reliable))
  expect_true(any(is.na(est$trt[!est$reliable])))
  expect_true(all(is.finite(est$trt[est$reliable])))
})

test_that("chestr_point returns skipped=TRUE below threshold", {
  skip_if_not_installed("survival")
  set.seed(7)
  n <- 40
  dat <- data.frame(
    st = rexp(n, 0.2),
    trt = rbinom(n, 1, 0.5),
    b1 = rnorm(n)
  )
  dat$cen <- as.numeric(dat$st < 2)
  dat$st[dat$cen == 0] <- 2
  base <- survival::coxph(survival::Surv(st, cen) ~ trt, data = dat)
  res <- suppressWarnings(
    chestr_point(base, dat["b1"], x = max(dat$b1) + 5,
                 method = "legacy", kern.adj = 10, min_events_per_df = 10)
  )
  expect_true(res$skipped)
  expect_null(res$fit)
})

test_that("plot.chestr uses stored base and biom", {
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("scales")
  set.seed(11)
  n <- 40
  dat <- data.frame(st = rexp(n, 0.2), trt = rbinom(n, 1, 0.5), b1 = rnorm(n))
  dat$cen <- 1L
  base <- survival::coxph(survival::Surv(st, cen) ~ trt, data = dat)
  cr <- chestr(base, dat$b1, grid.size = 3, min_events_per_df = 0)
  expect_s3_class(cr, "chestr")
  expect_false(inherits(cr, "data.frame"))
  p <- plot(cr, trt.param = "trt", pts = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot.chestr 2d returns ggplot", {
  skip_if_not_installed("survival")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("scales")
  set.seed(12)
  n <- 50
  dat <- data.frame(
    st = rexp(n, 0.2), trt = rbinom(n, 1, 0.5),
    b1 = rnorm(n), b2 = rnorm(n)
  )
  dat$cen <- 1L
  base <- survival::coxph(survival::Surv(st, cen) ~ trt, data = dat)
  cr <- chestr(base, dat[, c("b1", "b2")], grid.size = 4, min_events_per_df = 0)
  p <- plot(cr, trt.param = "trt", col.scale = c(-2, 2))
  expect_s3_class(p, "ggplot")
})

test_that("as.data.frame.chestr returns estimates", {
  skip_if_not_installed("survival")
  set.seed(3)
  n <- 50
  dat <- data.frame(
    st = rexp(n, 0.15),
    trt = rbinom(n, 1, 0.5),
    b1 = rnorm(n)
  )
  dat$cen <- as.numeric(dat$st < 3)
  dat$st[dat$cen == 0] <- 3
  base <- survival::coxph(survival::Surv(st, cen) ~ trt, data = dat)
  cr <- chestr(base, dat$b1, grid.size = 5, method = "distance",
               min_events_per_df = 0)
  df <- as.data.frame(cr)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), nrow(cr$estimates))
  expect_true(all(c("ess_events", "events_per_df", "reliable") %in% names(df)))
  expect_true(all(df$reliable))
})

test_that("chestr_statistics and chestr_test run", {
  skip_if_not_installed("survival")
  set.seed(21)
  n <- 60
  dat <- data.frame(
    st = rexp(n, 0.2),
    trt = rbinom(n, 1, 0.5),
    b1 = rnorm(n),
    b2 = rnorm(n)
  )
  dat$cen <- as.numeric(dat$st < 3)
  dat$st[dat$cen == 0] <- 3
  base <- survival::coxph(survival::Surv(st, cen) ~ trt, data = dat)
  cr <- chestr(base, dat[, c("b1", "b2")], grid.size = 4,
               min_events_per_df = 0, treat_term = "trt")
  expect_equal(cr$treat_term, "trt")
  expect_true(!is.null(cr$data))
  st <- chestr_statistics(cr)
  expect_true(is.finite(st$T_L2))
  expect_true(is.finite(st$T_MAX))
  expect_equal(st$treat_term, "trt")

  tst <- suppressWarnings(
    chestr_test(cr, B = 9, seed = 21, reliable_only = FALSE)
  )
  expect_s3_class(tst, "chestr_test")
  expect_equal(tst$treat_var, "trt")
  expect_true(tst$p_L2 >= 1 / (1 + tst$n_perm_used))
  expect_true(tst$p_L2 <= 1)
  expect_true(tst$p_MAX >= 1 / (1 + tst$n_perm_used))
  expect_equal(length(tst$null_T_L2), tst$n_perm_used)
})

test_that("chestr_test infers factor treat_var", {
  skip_if_not_installed("survival")
  set.seed(22)
  n <- 50
  dat <- data.frame(
    st = rexp(n, 0.2),
    treat_f = factor(sample(c("A", "B"), n, replace = TRUE)),
    b1 = rnorm(n)
  )
  dat$cen <- 1L
  base <- survival::coxph(survival::Surv(st, cen) ~ treat_f, data = dat)
  term <- grep("^treat_f", names(coef(base)), value = TRUE)[1]
  cr <- chestr(base, dat$b1, grid.size = 4, min_events_per_df = 0,
               treat_term = term)
  tst <- suppressWarnings(
    chestr_test(cr, B = 5, seed = 22, reliable_only = FALSE)
  )
  expect_equal(tst$treat_var, "treat_f")
  expect_equal(tst$treat_term, term)
})
