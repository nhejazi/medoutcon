# medoutcon 0.2.1

* Fixes issue affecting weighted TMLEs that introduced itself during
  previous update to `est_tml()`.

# medoutcon 0.2.0

* Added support for a semiparametric correction for outcome-dependent two-phase
  sampling designs with known or estimated sampling weights.
* Tightened sanity checks for estimation of natural direct and indirect effects
  by requiring that EIF scores related to intermediate confounders uniformly be
  zero when `Z = NULL` is specified.
* The above required minor changes to `est_tml()` so as to avoid fluctuation
  models for the TMLE step from updating the intermediate confounder nuisance
  components (i.e., `q_prime_Z_one`, `q_prime_Z_natural`) and causing numerical
  issues that violate the above internal checks.

# medoutcon 0.1.5

* For user clarity, the name of the argument for providing externally computed
  observation-level weights has changed (from `ext_weights`) to `svy_weights`.
* Support for the natural direct and indirect effects has been added, requiring
  the addition of the new internal argument `effect_type` across functions for
  estimation, including `cv_eif()`, `est_onestep()`, and `est_tml()`. When
  `Z = NULL` is set in `medoutcon()`, a natural effect estimate corresponding to
  the argument `effect` is returned instead of an interventional effect.
* The `summary()` and `print()` methods have been updated to allow handling of
  natural effects and counterfactual means under arbitrary contrasts.

# medoutcon 0.1.0

* An initial public release of this package, version 0.1.0, which includes
  support for external observation-level weights.
