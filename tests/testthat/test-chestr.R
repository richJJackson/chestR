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
  cr <- chestr(base, dat[, c("b1", "b2")], grid.size = 5,
               method = "legacy", kern.adj = 2)
  expect_true("trt" %in% names(cr))
  expect_true("eff.event" %in% names(cr))
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
  cr <- chestr(base, b, grid.size = 3, method = "legacy")
  expect_true(all(c("biom1", "biom2") %in% names(cr)))
  expect_true("x" %in% names(cr))
})
