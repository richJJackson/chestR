# chestR 0.1.0

* Initial packaged release of kernel-weighted Cox regression (`chestr`) for
  exploring treatment effect heterogeneity on a biomarker grid.
* Implements Kish effective-sample-size (ESS) events-per-df safeguard
  (`min_events_per_df`, default 10). Grid points below the threshold are
  skipped with `NA` coefficients and flagged via `reliable`.
* Exported API: `chestr()`, `chestr_point()`, `grid_biom()`, `plot_chestr()`
  (alias `plot.chestr()`), `chestr_statistics()`, `chestr_test()`.
* Internal helpers: `chestr_weights()`, `kish_ess()`, `euc.dist()`.
* Weighting methods: `"distance"` (default), `"square_distance"`, `"legacy"`.
* Global permutation test of constant treatment effect over biomarker space
  (`T_L2`, `T_MAX`) via `chestr_test()`.
* Vignette `chestr-workflow` demonstrates the global Cox → chestr → plot
  workflow (and optional permutation test) on simulated survival data.
