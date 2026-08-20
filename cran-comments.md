## Test environments

* Local: macOS (darwin), R 4.4.2 — `testthat::test_local()` and
  `R CMD check --no-manual --as-cran` during package prep.
* Not yet: win-builder, R-hub, or full multi-OS CRAN checks.

## R CMD check results

On the prep machine: **0 ERRORS, 0 WARNINGS**, with expected NOTES only:

* CRAN incoming feasibility: New submission.
* Future file timestamps: unable to verify current time (offline/clock check).

Vignette `chestr-workflow` builds from simulated data only.
A full CRAN-style check must still be re-run on a release tarball before
submission (see `CRAN_CHECKLIST.md`).

## Downstream dependencies

There are currently no reverse dependencies on CRAN.

## Notes for CRAN

* First submission of **chestR**.
* Maintainer will confirm email promptly after submission.
