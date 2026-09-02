# chestR

Kernel-weighted Cox regression for exploring treatment effect heterogeneity and
candidate predictive biomarkers.

## Installation

```r
# install.packages("devtools")
devtools::install_github("richJJackson/chestR")
```

Or install from a local checkout:

```r
devtools::install("path/to/chestR")
```

## Quick start

```r
library(survival)
library(chestR)

# Fit a global Cox model
base <- coxph(Surv(time, status) ~ treatment + covariate, data = mydata)

# Local estimates over a biomarker grid
cr <- chestr(base, mydata[, c("biom1", "biom2")], grid.size = 25,
             treat_term = "treatment")

# Visualise local treatment effect
plot(cr, trt.param = "treatment")

# Optional permutation test of constant treatment effect
# tst <- chestr_test(cr, B = 99, seed = 1)
```

See `vignette("chestr-workflow", package = "chestR")` after install,
or `inst/examples/simulation.R` for a longer simulation script.

## Development

Open `chestR.Rproj` in RStudio, then:

```r
devtools::load_all()
devtools::test()
devtools::document()
devtools::check()
```

## Related work

- Bonetti M, Gelber RD (2000). *Statistics in Medicine*.
- Liu Y, Lu W, Chen G (2015). Local partial-likelihood test. *Statistics in Medicine*.

## License

MIT
