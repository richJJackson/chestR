# CRAN submission checklist (chestR)

Do **not** submit until the items below are complete. This file is a prep
checklist only.

## Before `R CMD check --as-cran`

- [ ] Confirm maintainer email (`richj23@liverpool.ac.uk`) is monitored and
      will receive CRAN confirmation / follow-up mail.
- [ ] Run a full local check: `R CMD check --as-cran chestR_*.tar.gz`
      (or `devtools::check(cran = TRUE)`).
- [ ] Build and inspect vignettes: `devtools::build_vignettes()` / check that
      `chestr-workflow` builds cleanly without external data.
- [ ] Spell-check Description / Title / Rd pages (`devtools::spell_check()`).
- [ ] Confirm `LICENSE` / `LICENSE.md` handling matches MIT + file LICENSE.
- [ ] Ensure no non-ASCII issues, no unneeded files in the tarball
      (`.Rbuildignore` covers `.Rproj`, git, etc.).
- [ ] Decide whether `NEWS.md` should be installed (`inst/NEWS.md` or leave as
      package root; CRAN accepts root `NEWS.md`).
- [ ] Review examples: prefer `\dontrun{}` only where needed; vignette carries
      the runnable demo.
- [ ] Check reverse dependencies if any appear later (`revdepcheck`).

## Submission steps (when ready)

1. Bump `Version` if needed and update `NEWS.md`.
2. `devtools::build()` → `R CMD check --as-cran` on the tarball.
3. `devtools::submit_cran()` (or upload via CRAN web form).
4. Reply to the automated maintainer confirmation email promptly.
5. Watch for CRAN comments and fix promptly (usually within days).

## Known remaining items for first CRAN submission

- Maintainer email confirmation at submit time.
- Full `--as-cran` check on a clean machine / CRAN-like environment
  (Windows + macOS binaries after acceptance).
- Optional: ORCID / funding acknowledgements if desired in `DESCRIPTION`.
- Optional: CITATION file once a companion paper is available.
