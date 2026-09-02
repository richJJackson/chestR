## Test environments

* Local: macOS (darwin / aarch64), R 4.4.2 —
  `devtools::check(cran = TRUE)` and `R CMD check --as-cran` on
  `chestR_0.1.0.tar.gz` (2026-08-21).
* Not yet: win-builder, R-hub, or full multi-OS CRAN checks.

## R CMD check results

**0 ERRORS, 0 WARNINGS, 2 NOTES**

NOTES:

* `checking CRAN incoming feasibility ... NOTE`
  New submission (expected for a first CRAN upload).

* `checking for future file timestamps ... NOTE`
  `unable to verify current time`
  (common local/network clock-server check; not a package defect)

Vignette `chestr-workflow` builds from simulated data only.
Unit tests: 47 passed.
Spell check: clean (`Language: en-GB`, `inst/WORDLIST`).

## Downstream dependencies

There are currently no reverse dependencies on CRAN.

## Notes for CRAN

* First submission of **chestR**.
* Maintainer will confirm email promptly after submission.
