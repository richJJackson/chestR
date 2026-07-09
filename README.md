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
cr <- chestr(base, mydata$biomarker, grid.size = 25)

# Visualise local treatment effect (2 biomarkers)
plot_chestr(cr, trt.param = "treatment", base = base, biom = mydata[, c("biom1", "biom2")])

# After library(chestR), plot.chestr() is also available as a backward-compatible alias
```

See `inst/examples/simulation.R` for a self-contained simulation example.

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
