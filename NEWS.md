# MoCo 0.1.0

* Initial CRAN release.
* `moco()`: one-step and cross-fitted estimation of motion-controlled group
  means and contrasts, with HAL, log-normal GLM, or generalized-gamma GAMLSS
  conditional motion densities and Super Learner or GLM nuisance regressions.
* `support_rule` sets where the tolerable-motion density
  p(m | Delta_M = 1, A = 0, X) is zero: `"delta_m"` (default; every
  `Delta_M = 0` participant for GLM and GAMLSS) or `"observed_range"` (only
  above the largest motion among `Delta_M = 1` participants), for motion
  inclusion rules that are not a cutoff on motion.
* Motion densities used as denominators are bounded below by an absolute
  `density_floor` (default `sqrt(.Machine$double.eps)`); an optional
  `density_ratio_cap` (default `Inf`) is available for sensitivity analyses.
  Non-finite efficient-influence-function values stop the fit.
* GAMLSS densities use penalized splines for named continuous covariates, a
  configurable convergence criterion (`gamlss_c_crit`), and an automatic
  fallback to the RS, mixed and CG algorithms when the requested optimizer
  fails or does not converge; the algorithm used is recorded in the output.
* Repeated cross-fitting seeds (`seed_rgn`) are combined by averaging
  estimates and participant-level efficient influence functions before the
  covariance and simultaneous critical values are computed.
* `hypo_test()`: EIF-based simultaneous confidence-band tests controlling the
  family-wise error rate across outcomes.
* `plot_moco()`: optional cortical-surface plots of parcel-wise results
  (requires `ciftiTools` and Connectome Workbench).
* Analyses without motion data (`M = NULL`, optionally `X = NULL`) follow
  Nebel et al. (2022).
* Simulated example data (`data`), examples, tests and a citation.
