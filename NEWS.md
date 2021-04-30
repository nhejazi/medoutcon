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
