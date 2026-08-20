test_that("grid_biom returns expected dimensions", {
  g <- grid_biom(data.frame(x = 1:5, y = 11:15), grid.size = 4)
  expect_equal(nrow(g), 16)
  expect_equal(names(g), c("x", "y"))
})

test_that("chestr returns coefficients and effective events", {
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
  cr <- suppressWarnings(
    chestr(base, dat[, c("b1", "b2")], grid.size = 5,
           method = "legacy", kern.adj = 2, min_events_per_df = 5)
  )
  expect_true("trt" %in% names(cr))
  expect_true("eff.event" %in% names(cr))
  expect_true("ess_events" %in% names(cr))
  expect_equal(nrow(cr), 25)
  expect_true(any(grepl("\\.se$", names(cr))))
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
  expect_true(all(c("biom1", "biom2") %in% names(cr)))
  expect_true("x" %in% names(cr))
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
  expect_true("reliable" %in% names(cr))
  expect_true("ess_events" %in% names(cr))
  expect_true(any(cr$reliable))
  expect_true(any(!cr$reliable))
  expect_true(any(is.na(cr$trt[!cr$reliable])))
  expect_true(all(is.finite(cr$trt[cr$reliable])))
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


test_that("plot.chestr is available as plot() method", {
  expect_true(exists("plot.chestr", mode = "function"))
  expect_equal(names(formals(plot.chestr)), c("x", "..."))
})

test_that("chestr output has S3 class chestr", {
  skip_if_not_installed("survival")
  set.seed(11)
  n <- 40
  dat <- data.frame(st = rexp(n, 0.2), trt = rbinom(n, 1, 0.5), b1 = rnorm(n))
  dat$cen <- 1L
  base <- survival::coxph(survival::Surv(st, cen) ~ trt, data = dat)
  cr <- chestr(base, dat$b1, grid.size = 3, min_events_per_df = 0)
  expect_s3_class(cr, "chestr")
  expect_s3_class(cr, "data.frame")
})

test_that("distance method returns reliable ESS columns", {
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
  expect_true(all(c("ess_events", "events_per_df", "reliable") %in% names(cr)))
  expect_true(all(cr$reliable))
})
