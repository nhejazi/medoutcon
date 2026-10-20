# Efficient estimation of natural and interventional (in)direct effects

Efficient estimation of natural and interventional (in)direct effects

## Usage

``` r
medoutcon(
  W,
  A,
  Z,
  M,
  Y,
  R = rep(1, length(Y)),
  obs_weights = rep(1, length(Y)),
  svy_weights = NULL,
  two_phase_weights = rep(1, length(Y)),
  effect = c("direct", "indirect", "pm"),
  contrast = NULL,
  g_learners = sl3::Lrnr_glm_fast$new(),
  h_learners = sl3::Lrnr_glm_fast$new(),
  b_learners = sl3::Lrnr_glm_fast$new(),
  q_learners = sl3::Lrnr_glm_fast$new(),
  r_learners = sl3::Lrnr_glm_fast$new(),
  u_learners = sl3::Lrnr_hal9001$new(),
  v_learners = sl3::Lrnr_hal9001$new(),
  d_learners = sl3::Lrnr_glm_fast$new(),
  estimator = c("tmle", "onestep"),
  estimator_args = list(cv_folds = 10L, cv_strat = FALSE, strat_pmin = 0.1, max_iter =
    10L, tiltmod_tol = 5),
  g_bounds = c(0.005, 0.995)
)
```

## Arguments

- W:

  A `matrix`, `data.frame`, or similar object corresponding to a set of
  baseline covariates.

- A:

  A `numeric` vector corresponding to a treatment variable. The
  parameter of interest is defined as a location shift of this quantity.

- Z:

  A `numeric` vector corresponding to an intermediate confounder
  affected by treatment (on the causal pathway between the intervention
  A, mediators M, and outcome Y, but unaffected itself by the
  mediators). When set to `NULL`, the natural (in)direct effects are
  estimated.

- M:

  A `numeric` vector, `matrix`, `data.frame`, or similar corresponding
  to a set of mediators (on the causal pathway between the intervention
  A and the outcome Y).

- Y:

  A `numeric` vector corresponding to an outcome variable.

- R:

  A `logical` vector indicating whether a sampled observation's mediator
  was measured via a two-phase sampling design. Defaults to a vector of
  ones, indicating that two-phase sampling was not performed.

- obs_weights:

  A `numeric` vector of observation-level weights. The default is to
  give all observations equal weighting.

- svy_weights:

  A `numeric` vector of observation-level weights that have been
  computed externally, such as survey sampling weights. Such weights are
  used in the construction of re-weighted efficient estimators.

- two_phase_weights:

  A `numeric` vector of known observation-level weights corresponding to
  the inverse probability of the mediator being measured. Defaults to a
  vector of ones.

- effect:

  A `character` indicating whether to compute the direct or the indirect
  effects of \<https://doi.org/10.1093/biomet/asaa085\>. This is ignored
  when the argument `contrast` is provided. By default, the direct
  effect is estimated.

- contrast:

  A `numeric` double indicating the two values of the intervention `A`
  to be compared. The default value of `NULL` has no effect, as the
  value of the argument `effect` is instead used to define the
  contrasts. To override `effect`, provide a `numeric` double vector,
  giving the values of a' and a\*, e.g., `c(0, 1)`.

- g_learners:

  A [`Stack`](https://tlverse.org/sl3/reference/Stack.html) object, or
  other learner class (inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html)),
  containing instantiated learners from sl3; used to fit a model for the
  propensity score.

- h_learners:

  A [`Stack`](https://tlverse.org/sl3/reference/Stack.html) object, or
  other learner class (inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html)),
  containing instantiated learners from sl3; used to fit a model for a
  parameterization of the propensity score that conditions on the
  mediators.

- b_learners:

  A [`Stack`](https://tlverse.org/sl3/reference/Stack.html) object, or
  other learner class (inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html)),
  containing instantiated learners from sl3; used to fit a model for the
  outcome regression.

- q_learners:

  A [`Stack`](https://tlverse.org/sl3/reference/Stack.html) object, or
  other learner class (inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html)),
  containing instantiated learners from sl3; used to fit a model for a
  nuisance regression of the intermediate confounder, conditioning on
  the treatment and potential baseline covariates.

- r_learners:

  A [`Stack`](https://tlverse.org/sl3/reference/Stack.html) object, or
  other learner class (inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html)),
  containing instantiated learners from sl3; used to fit a model for a
  nuisance regression of the intermediate confounder, conditioning on
  the mediators, the treatment, and potential baseline confounders.

- u_learners:

  A [`Stack`](https://tlverse.org/sl3/reference/Stack.html) object, or
  other learner class (inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html)),
  containing instantiated learners from sl3; used to fit a
  pseudo-outcome regression required for in the efficient influence
  function.

- v_learners:

  A [`Stack`](https://tlverse.org/sl3/reference/Stack.html) object, or
  other learner class (inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html)),
  containing instantiated learners from sl3; used to fit a
  pseudo-outcome regression required for in the efficient influence
  function.

- d_learners:

  A [`Stack`](https://tlverse.org/sl3/reference/Stack.html) object, or
  other learner class (inheriting from
  [`Lrnr_base`](https://tlverse.org/sl3/reference/Lrnr_base.html)),
  containing instantiated learners from sl3; used to fit an initial
  efficient influence function regression when computing the efficient
  influence function in a two-phase sampling design.

- estimator:

  The desired estimator of the direct or indirect effect (or
  contrast-specific parameter) to be computed. Both an efficient
  one-step estimator using cross-fitting and a cross-validated targeted
  minimum loss estimator (TMLE) are available. The default is the TML
  estimator.

- estimator_args:

  A `list` of additional arguments passed (via `...`) to the function
  call for the specified estimator. The default is chosen so as to allow
  the number of folds used to compute the one-step or TML estimators to
  be adjusted and for stratified cross-validation to be used in cases of
  rare outcomes. In the case of the TML estimator, the number of update
  (fluctuation) iterations is limited, and a tolerance is included for
  the updates introduced by tilting (fluctuation) models.

- g_bounds:

  A `numeric` vector containing two values, the first being the minimum
  allowable estimated propensity score value and the second being the
  maximum allowable for estimated propensity scores.

## Examples

``` r
# here, we show one-step and TML estimates of the interventional direct
# effect; the indirect effect can be evaluated by a straightforward change
# to the penultimate argument. the natural direct and indirect effects can
# be evaluated by omitting the argument Z (inappropriate in this example).
# create data: covariates W, exposure A, post-exposure-confounder Z,
#              mediator M, outcome Y
n_obs <- 200
w_1 <- rbinom(n_obs, 1, prob = 0.6)
w_2 <- rbinom(n_obs, 1, prob = 0.3)
w <- as.data.frame(cbind(w_1, w_2))
a <- as.numeric(rbinom(n_obs, 1, plogis(rowSums(w) - 2)))
z <- rbinom(n_obs, 1, plogis(rowMeans(-log(2) + w - a) + 0.2))
m_1 <- rbinom(n_obs, 1, plogis(rowSums(log(3) * w + a - z)))
m_2 <- rbinom(n_obs, 1, plogis(rowSums(w - a - z)))
m <- as.data.frame(cbind(m_1, m_2))
y <- rbinom(n_obs, 1, plogis(1 / (rowSums(w) - z + a + rowSums(m))))

# one-step estimate of the interventional direct effect
os_de <- medoutcon(
  W = w, A = a, Z = z, M = m, Y = y,
  effect = "direct",
  estimator = "onestep"
)

# TML estimate of the interventional direct effect
# NOTE: improved variance estimate and de-biasing from targeting procedure
tmle_de <- medoutcon(
  W = w, A = a, Z = z, M = m, Y = y,
  effect = "direct",
  estimator = "tmle"
)
```
