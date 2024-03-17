# medoutcon 0.2.3

* Added a new named argument `cv_stratify` to `est_onestep()` and `est_tml()`
  and to the `estimator_args` list-argument in `medoutcon()`, which allows for
  stratified folds to be generated for cross-fitting (by passing these to the
  `strata_ids` argument of `make_folds()` from the `origami` package). This is
  also triggered by an override in `est_onestep()` and `est_tml()` when the
  proportion of detected cases is less than 0.1, a heuristic for rare outcomes.
* Increased the default number of folds for cross-fitting from 5 to 10, setting
  `cv_folds = 10L` in named arguments to `est_onestep()` and `est_tml()` and to
  the `estimator_args` list-argument in `medoutcon()`.
* Changed default propensity score truncation bounds specified in `g_bounds` to
  `c(0.005, 0.995)` from `c(0.01, 0.99)` (in v0.22), based on sanity checks and 
  manual experimentation.
* Wrapped instances of `sl3_Task()` in which `outcome_type = "continuous"` is
  specified in `suppressWarnings()` to sink warnings when the outcome variable
  for a given nuisance estimation task fails `sl3`'s check for continuous-ness.

# medoutcon 0.2.2

* Change iterative targeting procedures in `est_tml()` to use `glm2::glm2`
  instead of `stats::glm` to avoid issues with erratic  IRLS in the latter;
  see <https://journal.r-project.org/archive/2011/RJ-2011-012/> for details.
* Changed default propensity score truncation bounds specified in `g_bounds` by
  an order of magnitude, from `c(0.001, 0.999)` to `c(0.01, 0.99)`, to mitigate
  potential stability issues.

# medoutcon 0.2.1

* Fixes bug in weighted TMLEs introduced during prior update to `est_tml()`.

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
