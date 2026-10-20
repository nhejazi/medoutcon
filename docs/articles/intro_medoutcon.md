# Efficient causal mediation analysis with the natural and interventional effects

## Background and Motivations

An exposure of interest often affects an outcome directly, or indirectly
by the mediation of some intermediate variables. Identifying and
quantifying the mechanisms underlying causal effects is an increasingly
popular endeavor in public health, medicine, and the social sciences, as
knowledge of such mechanisms can improve understanding of both *why and
how* treatments can be effective. Such mechanistic knowledge may be
arguably even more important in cases where treatments result in
unanticipated ineffective or even harmful effects.

Traditional techniques for mediation analysis fare poorly in the face of
intermediate confounding. Classical parameters like the natural
(in)direct effects face a lack of identifiability in cases where
mediator-outcome (i.e., intermediate) confounders affected by exposure
complicate the relationship between the exposure, mediators, and
outcome. Dı́az et al. (2020) provide a theoretical and computational
study of the properties of newly developed interventional (in)direct
effect estimands within the non-parametric statistical model. Among
their contributions, Dı́az et al. (2020)

- derive the efficient influence function (EIF), an key object in
  semiparametric efficiency theory;
- use the EIF to develop two asymptotically optimal, non-parametric
  estimators, each of which is capable of leveraging machine learning
  for the estimation of nuisance parameters; and
- present theoretical conditions under which their proposed estimators
  are consistent, multiply robust, and efficient.

## Problem Setup and Notation

The problem addressed by the work of Dı́az et al. (2020) may be
represented by the following nonparametric structural equation model
(NPSEM):
``` math
\begin{align*}
  W &= f_W(U_W); A = f_A(W, U_A); Z=f_Z(W, A, U_Z);\\ \nonumber
  M &= f_M(W, A, Z, U_M); Y = f_Y(W, A, Z, M, U_Y).
\end{align*}
```
In the NPSEM, $`W`$ denotes a vector of observed pre-treatment
covariates, $`A`$ denotes a categorical treatment variable, $`Z`$
denotes an intermediate confounder affected by treatment, $`M`$ denotes
a (possibly multivariate) mediator, and $`Y`$ denotes a continuous or
binary outcome. The vector of exogenous factors
$`U=(U_W,U_A,U_Z,U_M,U_Y)`$, and the functions $`f`$, are assumed
deterministic but unknown. Importantly, the NPSEM encodes a
time-ordering between these variables and allows the evaluation of
counterfactual quantities defined by intervening on a set of nodes of
the NPSEM. The observed data unit can be represented by the random
variable $`O = (W, A, Z, M, Y)`$; we consider access to
$`O_1, \ldots, O_n`$, a sample of $`n`$ i.i.d. observations of $`O`$.

Dı́az et al. (2020) additionally define the following parameterizations,
familiarity with which will be useful for using the [`medoutcon` `R`
package](https://github.com/nhejazi/medoutcon). In particular, these
authors define $`g(a \mid w)`$ as the probability mass function of
$`A = a`$ conditional on $`W = w`$ and use $`h(a \mid m, w)`$ to denote
the probability mass function of $`A
= a`$ conditional on $`(M, W) = (m, w)`$. Further, Dı́az et al. (2020)
use $`b(a, z, m, w)`$ to denote the outcome regression function
$`\mathbb{E}(Y \mid A
= a, Z = z, M = m, W = w)`$, as well as $`q(z \mid a,w)`$ and
$`r(z \mid a, m, w)`$ to denote the corresponding conditional densities
of $`Z`$.

## Interventional (In)Direct Effects

Dı́az et al. (2020) define the *total effect* of $`A`$ on $`Y`$ in terms
of a contrast between two user-supplied values
$`a', a^{\star} \in \mathcal{A}`$. Examination of the NPSEM reveals that
there are four paths involved in this effect, namely
$`A \rightarrow Y`$, $`A \rightarrow M \rightarrow Y`$, $`A
\rightarrow Z \rightarrow Y`$, and
$`A \rightarrow Z \rightarrow M \rightarrow Y`$. Mediation analysis has
classically considered the *natural direct effect* (NDE) and the
*natural indirect effect* (NIE), which are defined as
$`\mathbb{E}_c(Y_{a', M_{a^{\star}}} - Y_{a^{\star}, M_{a^{\star}}})`$
and $`\mathbb{E}_c(Y_{a',M_{a'}} - Y_{a',M_{a^{\star}}})`$,
respectively. The natural direct effect measures the effect through
paths *not* involving the mediator ($`A \rightarrow Y`$ and
$`A \rightarrow Z \rightarrow Y`$), whereas the natural indirect effect
measures the effect through paths involving the mediator
($`A \rightarrow M \rightarrow Y`$ and
$`A \rightarrow Z \rightarrow M \rightarrow
Y`$). As the sum of the natural direct and indirect effects equals the
average treatment effect $`\mathbb{E}_c(Y_1-Y_0)`$, this effect
decomposition is appealing. Unfortunately, the natural direct and
indirect effects are not generally identified in the presence of an
intermediate confounder affected by treatment.

To circumvent this issue, Dı́az et al. (2020) define the direct and
indirect effects using stochastic interventions on the mediator,
following a strategy previously outlined by VanderWeele et al. (2014)
and Rudolph et al. (2017), among others. Let $`G_a`$ denote a random
draw from the conditional distribution of $`M_a`$ conditional on $`W`$.
Consider the effect of $`A`$ on $`Y`$ defined as the difference in
expected outcome in hypothetical worlds in which $`(A,M) = (a',
G_{a'})`$ versus $`(A,M) = (a^{\star}, G_{a^{\star}})`$ with probability
one, which may be decomposed into direct and indirect effects as follows
``` math
\begin{equation*}
\mathbb{E}_c(Y_{a', G_{a'}} - Y_{a^{\star}, G_{a^{\star}}}) =
  \underbrace{\mathbb{E}_c(Y_{a', G_{a'}} - Y_{a',
    G_{a^{\star}}})}_{\text{Indirect effect (through $M$)}} +
  \underbrace{\mathbb{E}_c(Y_{a', G_{a^{\star}}} - Y_{a^{\star},
      G_{a^{\star}}})}_{\text{Direct effect (not through $M$)}}.
\end{equation*}
```
Like the natural direct effect, this interventional direct effect
measures the effects through paths not involving the mediator. Likewise,
the interventional indirect effect measures the effect through paths
involving the mediator. Note, however, that natural and interventional
mediation effects have different interpretations. That is, the
interventional indirect effect measures the effect of fixing the
exposure at $`a'`$ while setting the mediator to a random draw
$`G_{a^{\star}}`$ from those with exposure $`a'`$ versus a random draw
$`G_{a'}`$ from those with exposure $`a^{\star}`$, given covariates
$`W`$. As is clear from the effect decomposition, the term
$`\theta_c = \mathbb{E}_c(Y_{a', G_{a^{\star}}})`$ is required for
estimation of both the interventional direct and indirect effects; thus,
Dı́az et al. (2020) focus on estimation of this quantity. Importantly, it
has been shown that $`\theta_c`$ is identified by the statistical
functional
``` math
\begin{equation*}
  \theta = \int b(a', z, m, w) q(z \mid a', w) p(m \mid a^{\star}, w)
    p(w) d\nu(w,z,m)
\end{equation*}
```
under a set of standard identifiability conditions (VanderWeele et al.
2014), which are further reviewed in Dı́az et al. (2020).

## Efficient Estimation

Dı́az et al. (2020) define two efficient estimators of their
interventional (in)direct effects. These are based on the one-step
estimation and targeted minimum loss (TML) estimation frameworks,
respectively. Briefly, both estimation strategies proceed in two stages,
starting by first constructing initial estimates of the nuisance
parameters present in the EIF, then proceeding to apply distinct
bias-correction strategies in their second stages. Both estimation
strategies require an assumption about the behavior of initial
estimators of the nuisance parameters (specifically, that these lie in a
Donsker class); however, the need for such an assumption may be avoided
by making use of cross-validation in the fitting fo initial estimators.
The `medoutcon` `R` package requires the use of cross-validation in the
construction of these initial estimates, resulting in cross-fitted
one-step and and cross-validated TML estimators (Klaassen 1987; Zheng
and van der Laan 2011; Chernozhukov et al. 2018).

The one-step estimator $`\hat{\theta}_{\text{os}}`$ is constructed by
adding the empirical mean of the EIF (evaluated at initial estimates of
the nuisance parameters) to the substitution estimator. By constrast,
the TML estimator $`\hat{\theta}_{\text{tmle}}`$ updates the components
of the substitution estimator via logistic tilting models formulated to
ensure that relevant score equations appearing in the EIF are
(approximately) solved. While the estimators are asymptotically
equivalent, TML estimators have been shown to exhibit superior
finite-sample performance, making them potentially more reliable than
one-step estimators. For the exact form of the EIF as well as those of
the one-step and TML estimators, consult Dı́az et al. (2020).

## Data Analysis Example

### Setting up the data example

Now, we’ll take a look at how to estimate the interventional direct and
indirect effects using a simulated data example. Dı́az et al. (2020)
illustrate the use of their estimators of these effects in an
application in which they seek to elucidate the mechanisms behind the
unintended harmful effects that a housing intervention had on adolescent
girls’ risk behavior.

First, let’s load a few required packages and set a seed for our
simulation.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`data.table`](https://r-datatable.com)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`medoutcon`](https://github.com/nhejazi/medoutcon)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`sl3`](https://tlverse.org/sl3)`)`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``75681``)`\
`n_obs`` ``<-`` ``500`

Next, we’ll generate a very simple simulated dataset. The function
`make_example_data`, defined below, generates three binary baseline
covariates $`W = (W_1, W_2, W_3)`$, a binary exposure variable $`A`$, a
single binary mediateor $`M`$ an intermediate confounder $`Z`$ that
affects the mediator $`M`$ and is itself affected by the exposure $`A`$,
and, finally, a binary outcome $`Y`$ that is a function of
$`(W, A, Z, M)`$.

\
`# produces a simple data set based on ca causal model with mediation`\
`make_example_data`` ``<-`` ``function``(``n_obs`` ``=`` ``1000``)`` ``{`\
`  ``## baseline covariates`\
`  ``w_1`` ``<-`` `[`rbinom`](https://rdrr.io/r/stats/Binomial.html)`(``n_obs``, ``1``, prob ``=`` ``0.6``)`\
`  ``w_2`` ``<-`` `[`rbinom`](https://rdrr.io/r/stats/Binomial.html)`(``n_obs``, ``1``, prob ``=`` ``0.3``)`\
`  ``w_3`` ``<-`` `[`rbinom`](https://rdrr.io/r/stats/Binomial.html)`(``n_obs``, ``1``, prob ``=`` `[`pmin`](https://rdrr.io/r/base/Extremes.html)`(``0.2`` ``+`` ``(``w_1`` ``+`` ``w_2``)`` ``/`` ``3``, ``1``)``)`\
`  ``w`` ``<-`` `[`cbind`](https://rdrr.io/r/base/cbind.html)`(``w_1``, ``w_2``, ``w_3``)`\
`  ``w_names`` ``<-`` `[`paste`](https://rdrr.io/r/base/paste.html)`(``"W"``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`ncol`](https://rdrr.io/r/base/nrow.html)`(``w``)``)``, sep ``=`` ``"_"``)`\
\
`  ``## exposure`\
`  ``a`` ``<-`` `[`as.numeric`](https://rdrr.io/r/base/numeric.html)`(`[`rbinom`](https://rdrr.io/r/stats/Binomial.html)`(``n_obs``, ``1``, `[`plogis`](https://rdrr.io/r/stats/Logistic.html)`(`[`rowSums`](https://rdrr.io/r/base/colSums.html)`(``w``)`` ``-`` ``2``)``)``)`\
\
`  ``## mediator-outcome confounder affected by treatment`\
`  ``z`` ``<-`` `[`rbinom`](https://rdrr.io/r/stats/Binomial.html)`(``n_obs``, ``1``, `[`plogis`](https://rdrr.io/r/stats/Logistic.html)`(`[`rowMeans`](https://rdrr.io/r/base/colSums.html)`(``-`[`log`](https://rdrr.io/r/base/Log.html)`(``2``)`` ``+`` ``w`` ``-`` ``a``)`` ``+`` ``0.2``)``)`\
\
`  ``## mediator -- could be multivariate`\
`  ``m`` ``<-`` `[`rbinom`](https://rdrr.io/r/stats/Binomial.html)`(``n_obs``, ``1``, `[`plogis`](https://rdrr.io/r/stats/Logistic.html)`(`[`rowSums`](https://rdrr.io/r/base/colSums.html)`(`[`log`](https://rdrr.io/r/base/Log.html)`(``3``)`` ``*`` ``w``[``, ``-``3``]`` ``+`` ``a`` ``-`` ``z``)``)``)`\
`  ``m_names`` ``<-`` ``"M"`\
\
`  ``## outcome`\
`  ``y`` ``<-`` `[`rbinom`](https://rdrr.io/r/stats/Binomial.html)`(``n_obs``, ``1``, `[`plogis`](https://rdrr.io/r/stats/Logistic.html)`(``1`` ``/`` ``(`[`rowSums`](https://rdrr.io/r/base/colSums.html)`(``w``)`` ``-`` ``z`` ``+`` ``a`` ``+`` ``m``)``)``)`\
\
`  ``## construct output`\
`  ``dat`` ``<-`` `[`as.data.table`](https://rdrr.io/pkg/data.table/man/as.data.table.html)`(`[`cbind`](https://rdrr.io/r/base/cbind.html)`(``w ``=`` ``w``, a ``=`` ``a``, z ``=`` ``z``, m ``=`` ``m``, y ``=`` ``y``)``)`\
`  `[`setnames`](https://rdrr.io/pkg/data.table/man/setattr.html)`(``dat``, `[`c`](https://rdrr.io/r/base/c.html)`(``w_names``, ``"A"``, ``"Z"``, ``m_names``, ``"Y"``)``)`\
`  `[`return`](https://rdrr.io/r/base/function.html)`(``dat``)`\
`}`\
\
`# set seed and simulate example data`\
`example_data`` ``<-`` ``make_example_data``(``n_obs``)`\
`w_names`` ``<-`` ``stringr``::`[`str_subset`](https://stringr.tidyverse.org/reference/str_subset.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``example_data``)``, ``"W"``)`\
`m_names`` ``<-`` ``stringr``::`[`str_subset`](https://stringr.tidyverse.org/reference/str_subset.html)`(`[`colnames`](https://rdrr.io/r/base/colnames.html)`(``example_data``)``, ``"M"``)`

Now, let’s take a quick look at our simulated data:

\
`# quick look at the data`\
[`head`](https://rdrr.io/r/utils/head.html)`(``example_data``)`

    ##      W_1   W_2   W_3     A     Z     M     Y
    ##    <num> <num> <num> <num> <num> <num> <num>
    ## 1:     0     0     1     0     0     0     1
    ## 2:     0     1     1     0     0     0     1
    ## 3:     1     1     0     0     1     1     1
    ## 4:     1     0     1     0     1     0     0
    ## 5:     0     0     0     1     1     1     1
    ## 6:     1     0     1     1     0     1     1

As noted above, all covariates in our dataset are binary; however, note
that this need not be the case for using our methodology — in
particular, the only current limitation is that the intermediate
confounder $`Z`$ must be binary when using our implemented TML estimator
of the (in)direct effects.

Using this dataset, we’ll proceed to estimate the interventional
(in)direct effects. In order to do so, we’ll need to estimate several
nuisance parameters, including the exposure mechanism $`g(A \mid W)`$, a
re-parameterized exposure mechanism that conditions on the mediators
$`h(A \mid M, W)`$, the outcome mechanism $`b(Y \mid M, Z, A, W)`$, and
two variants of the intermediate confounding mechanism
$`q(Z \mid A, W)`$ and $`r(Z \mid M, A, W)`$. In order to estimate each
of these nuisance parameters flexibly, we’ll rely on data adaptive
regression strategies in order to avoid the potential for (parametric)
model misspecification.

### Ensemble learning of nuisance functions

As we’d like to rely on flexible, data adaptive regression strategies
for estimating each of the nuisance parameters $`(g, h, b, q, r)`$, we
require a method for choosing among or combining the wide variety of
available regression strategies. For this, we recommend the use of the
Super Learner algorithm for ensemble machine learning (van der Laan et
al. 2007). The recently developed [`sl3` R
package](https://tlverse.org/sl3) (**coyle2020sl3?**) provides a unified
interface for deploying a wide variety of machine learning algorithms
(simply called *learners* in the `sl3` nomenclature) as well as for
constructing Super Learner ensemble models of such learners. For a
complete guide on using the `sl3` R package, consider consulting
<https://tlverse.org/sl3>, or <https://tlverse.org> (and
<https://github.com/tlverse>) for the `tlverse` ecosystem, of which
`sl3` is an integral part.

To construct an ensemble learner using a handful of popular machine
learning algorithms, we’ll first instantiate variants of learners from
the appropriate classes for each algorithm, and then create a Super
Learner ensemble via the `Lrnr_sl` class. Below, we demonstrate the
construction of an ensemble learner based on a modeling library
including an intercept model, a main-terms GLM, $`\ell_1`$-penalized
Lasso regression, an elastic net regression that equally weights the
$`\ell_1`$ and $`\ell_2`$ penalties, random forests (`ranger`), and the
highly adaptive lasso (HAL):

\
`# instantiate learners`\
`mean_lrnr`` ``<-`` `[`Lrnr_mean`](https://tlverse.org/sl3/reference/Lrnr_mean.html)`$``new``(``)`\
`fglm_lrnr`` ``<-`` `[`Lrnr_glm_fast`](https://tlverse.org/sl3/reference/Lrnr_glm_fast.html)`$``new``(``family ``=`` `[`binomial`](https://rdrr.io/r/stats/family.html)`(``)``)`\
`lasso_lrnr`` ``<-`` `[`Lrnr_glmnet`](https://tlverse.org/sl3/reference/Lrnr_glmnet.html)`$``new``(``alpha ``=`` ``1``, family ``=`` ``"binomial"``, nfolds ``=`` ``3``)`\
`enet_lrnr`` ``<-`` `[`Lrnr_glmnet`](https://tlverse.org/sl3/reference/Lrnr_glmnet.html)`$``new``(``alpha ``=`` ``0.5``, family ``=`` ``"binomial"``, nfolds ``=`` ``3``)`\
`rf_lrnr`` ``<-`` `[`Lrnr_ranger`](https://tlverse.org/sl3/reference/Lrnr_ranger.html)`$``new``(``num.trees ``=`` ``200``)`\
\
`# for HAL, use linear probability formulation, with bounding in unit interval`\
`hal_gaussian_lrnr`` ``<-`` `[`Lrnr_hal9001`](https://tlverse.org/sl3/reference/Lrnr_hal9001.html)`$``new``(`\
`  family ``=`` ``"gaussian"``,`\
`  fit_control ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`    max_degree ``=`` ``3``,`\
`    n_folds ``=`` ``3``,`\
`    use_min ``=`` ``TRUE``,`\
`    type.measure ``=`` ``"mse"`\
`  ``)`\
`)`\
`bound_lrnr`` ``<-`` `[`Lrnr_bound`](https://tlverse.org/sl3/reference/Lrnr_bound.html)`$``new``(``bound ``=`` ``1e-6``)`\
`hal_bounded_lrnr`` ``<-`` `[`Pipeline`](https://tlverse.org/sl3/reference/Pipeline.html)`$``new``(``hal_gaussian_lrnr``, ``bound_lrnr``)`\
\
`# create learner library and instantiate super learner ensemble`\
`lrnr_lib`` ``<-`` `[`Stack`](https://tlverse.org/sl3/reference/Stack.html)`$``new``(`\
`  ``mean_lrnr``, ``fglm_lrnr``, ``enet_lrnr``, ``lasso_lrnr``,`\
`  ``rf_lrnr``, ``hal_bounded_lrnr`\
`)`\
`sl_lrnr`` ``<-`` `[`Lrnr_sl`](https://tlverse.org/sl3/reference/Lrnr_sl.html)`$``new``(``learners ``=`` ``lrnr_lib``, metalearner ``=`` `[`Lrnr_nnls`](https://tlverse.org/sl3/reference/Lrnr_nnls.html)`$``new``(``)``)`

While we recommend the use of a Super Learner ensemble model like the
one constructed above in practice, such a library will be too
computationally intensive for our examples. To reduce computation time,
we construct a simpler library, using only a subset of the above
learning algorithms:

\
`# create simpler learner library and instantiate super learner ensemble`\
`lrnr_lib`` ``<-`` `[`Stack`](https://tlverse.org/sl3/reference/Stack.html)`$``new``(``mean_lrnr``, ``fglm_lrnr``, ``lasso_lrnr``, ``rf_lrnr``)`\
`sl_lrnr`` ``<-`` `[`Lrnr_sl`](https://tlverse.org/sl3/reference/Lrnr_sl.html)`$``new``(``learners ``=`` ``lrnr_lib``, metalearner ``=`` `[`Lrnr_nnls`](https://tlverse.org/sl3/reference/Lrnr_nnls.html)`$``new``(``)``)`

Having set up our ensemble learner, we’re now ready to estimate each of
the interventional effects using the efficient estimators exposed in the
`medoutcon` package.

### Estimating the direct effect

We’re now ready to estimate the interventional direct effect. This
direct effect is computed as a contrast between the interventions
$`(a' = 1, a^{\star} = 0)`$ and $`(a' = 0, a^{\star} = 0)`$. In
particular, our efficient estimators of the interventional direct effect
proceed by constructing estimators
$`\hat{\theta}(a' = 1, a^{\star} = 0)`$ and
$`\hat{\theta}(a' = 0, a^{\star} = 0)`$. Then, an efficient estimator of
the direct effect is available by application of the delta method, that
is, $`\hat{\theta}^{\text{DE}} =
\hat{\theta}(a' = 1, a^{\star} = 0) - \hat{\theta}(a' = 0, a^{\star} = 0)`$.
Applying the same principle to the EIF estimates, one can derive
variance estimates and construct asymptotically correct Wald-style
confidence intervals for $`\hat{\theta}^{\text{DE}}`$.

The `medoutcon` package makes the estimation task quite simple, as only
a single call to the eponymous `medoutcon` function is required. As
demonstrated below, we need only feed in each component of the observed
data $`O = (W, A, Z, M, Y)`$ (of which $`W`$ and $`M`$ can be
multivariate), specify the effect type, and the estimator. Additionally,
for each nuisance parameter we may specify a separate regression
function — in the examples below, we use the simpler Super Learner
ensemble constructed above for fitting each nuisance function, but this
need not be the case (i.e., different estimators may be used for each
nuisance function).

First, we examine the one-step estimator of the interventional direct
effect. Recall that the one-step estimator is constructed by adding the
mean of the EIF (evaluated at initial estimates of the nuisance
parameters) to the substitution estimator. As noted above, this is done
separately for each of the two contrasts $`(a' = 0, a^{\star} = 0)`$ and
$`(a' = 1, a^{\star} = 0)`$. Thus, the one-step estimator of this direct
effect is constructed by application of the delta method to each of the
one-step estimators (and EIFs) for these contrasts.

\
`# compute one-step estimate of the interventional direct effect`\
`os_de`` ``<-`` `[`medoutcon`](https://codex.nimahejazi.org/medoutcon/reference/medoutcon.md)`(`\
`  W ``=`` ``example_data``[``, ``..w_names``]``,`\
`  A ``=`` ``example_data``$``A``,`\
`  Z ``=`` ``example_data``$``Z``,`\
`  M ``=`` ``example_data``[``, ``..m_names``]``,`\
`  Y ``=`` ``example_data``$``Y``,`\
`  g_learners ``=`` ``sl_lrnr``,`\
`  h_learners ``=`` ``sl_lrnr``,`\
`  b_learners ``=`` ``sl_lrnr``,`\
`  q_learners ``=`` ``sl_lrnr``,`\
`  r_learners ``=`` ``sl_lrnr``,`\
`  effect ``=`` ``"direct"``,`\
`  estimator ``=`` ``"onestep"``,`\
`  estimator_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``cv_folds ``=`` ``2``)`\
`)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``os_de``)`

    ## # A tibble: 1 × 7
    ##    lwr_ci param_est upr_ci var_est eif_mean estimator param                
    ##     <dbl>     <dbl>  <dbl>   <dbl>    <dbl> <chr>     <chr>                
    ## 1 -0.0881    0.0369  0.162 0.00406 2.84e-17 onestep   direct_interventional

From the output of the summary method, we note that the one-step
estimate of the interventional direct effect
$`\hat{\theta}_{\text{os}}^{\text{DE}}`$ is 0.037, with 95% confidence
interval \[-0.088, 0.162\].

Next, let’s compare the one-step estimate to the TML estimate. Analogous
to the case of the one-step estimator, the TML estimator can be
evaluated via a single call to the `medoutcon` function:

\
`# compute targeted minimum loss estimate of the interventional direct effect`\
`tmle_de`` ``<-`` `[`medoutcon`](https://codex.nimahejazi.org/medoutcon/reference/medoutcon.md)`(`\
`  W ``=`` ``example_data``[``, ``..w_names``]``,`\
`  A ``=`` ``example_data``$``A``,`\
`  Z ``=`` ``example_data``$``Z``,`\
`  M ``=`` ``example_data``[``, ``..m_names``]``,`\
`  Y ``=`` ``example_data``$``Y``,`\
`  g_learners ``=`` ``sl_lrnr``,`\
`  h_learners ``=`` ``sl_lrnr``,`\
`  b_learners ``=`` ``sl_lrnr``,`\
`  q_learners ``=`` ``sl_lrnr``,`\
`  r_learners ``=`` ``sl_lrnr``,`\
`  effect ``=`` ``"direct"``,`\
`  estimator ``=`` ``"tmle"``,`\
`  estimator_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``cv_folds ``=`` ``2``, max_iter ``=`` ``5``)`\
`)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``tmle_de``)`

    ## # A tibble: 1 × 7
    ##    lwr_ci param_est upr_ci var_est eif_mean estimator param                
    ##     <dbl>     <dbl>  <dbl>   <dbl>    <dbl> <chr>     <chr>                
    ## 1 -0.0736    0.0453  0.164 0.00368  0.00237 tmle      direct_interventional

From the output of the summary method, we note that the TML estimate of
the interventional direct effect
$`\hat{\theta}_{\text{tmle}}^{\text{DE}}`$ is 0.045, with 95% confidence
interval \[-0.074, 0.164\]. Here, we recall that the TML estimator
generally exhibits better finite-sample performance than the one-step
estimator (van der Laan and Rose 2011, 2018), so the TML estimate is
likely to be more reliable in our modest sample size of $`n =`$ 500.

### Estimating the indirect effect

Estimation of the interventional indirect effect proceeds similarly to
the strategy discussed above for the corresponding direct effect. An
efficient estimator can be computed as a contrast between the
interventions $`(a' = 1,
a^{\star} = 0)`$ and $`(a' = 1, a^{\star} = 1)`$. Specifically, our
efficient estimators of the interventional indirect effect proceed by
constructing estimators $`\hat{\theta}(a' = 1, a^{\star} = 0)`$ and
$`\hat{\theta}(a' = 1,
a^{\star} = 1)`$. Then, application of the delta method yields an
efficient estimator of the indirect effect, that is,
$`\hat{\theta}^{\text{IE}} =
\hat{\theta}(a' = 1, a^{\star} = 0) - \hat{\theta}(a' = 1, a^{\star} = 1)`$.
The same principle may be applied to the EIF estimates to derive
variance estimates and construct asymptotically correct Wald-style
confidence intervals for $`\hat{\theta}^{\text{IE}}`$.

Now, we examine the one-step estimator of the interventional indirect
effect. The one-step estimator is constructed by adding the mean of the
EIF (evaluated at initial estimates of the nuisance parameters) to the
substitution estimator. As noted above, this is done separately for each
of the two contrasts $`(a' = 1, a^{\star} = 1)`$ and
$`(a' = 1, a^{\star} = 0)`$. Thus, the one-step estimator of this
indirect effect is constructed by application of the delta method to
each of the one-step estimators (and EIFs) for the contrasts.

\
`# compute one-step estimate of the interventional indirect effect`\
`os_ie`` ``<-`` `[`medoutcon`](https://codex.nimahejazi.org/medoutcon/reference/medoutcon.md)`(`\
`  W ``=`` ``example_data``[``, ``..w_names``]``,`\
`  A ``=`` ``example_data``$``A``,`\
`  Z ``=`` ``example_data``$``Z``,`\
`  M ``=`` ``example_data``[``, ``..m_names``]``,`\
`  Y ``=`` ``example_data``$``Y``,`\
`  g_learners ``=`` ``sl_lrnr``,`\
`  h_learners ``=`` ``sl_lrnr``,`\
`  b_learners ``=`` ``sl_lrnr``,`\
`  q_learners ``=`` ``sl_lrnr``,`\
`  r_learners ``=`` ``sl_lrnr``,`\
`  effect ``=`` ``"indirect"``,`\
`  estimator ``=`` ``"onestep"`\
`)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``os_ie``)`

    ## # A tibble: 1 × 7
    ##   lwr_ci param_est  upr_ci var_est  eif_mean estimator param                  
    ##    <dbl>     <dbl>   <dbl>   <dbl>     <dbl> <chr>     <chr>                  
    ## 1 -0.183    -0.112 -0.0421 0.00129 -4.80e-17 onestep   indirect_interventional

From the output of the summary method, we note that the one-step
estimate of the interventional indirect effect
$`\hat{\theta}_{\text{os}}^{\text{IE}}`$ is -0.112, with 95% confidence
interval \[-0.183, -0.042\].

As before, let’s compare the one-step estimate to the TML estimate.
Analogous to the case of the one-step estimator, the TML estimator can
be evaluated via a single call to the `medoutcon` function, as
demonstrated below

\
`# compute targeted minimum loss estimate of the interventional indirect effect`\
`tmle_ie`` ``<-`` `[`medoutcon`](https://codex.nimahejazi.org/medoutcon/reference/medoutcon.md)`(`\
`  W ``=`` ``example_data``[``, ``..w_names``]``,`\
`  A ``=`` ``example_data``$``A``,`\
`  Z ``=`` ``example_data``$``Z``,`\
`  M ``=`` ``example_data``[``, ``..m_names``]``,`\
`  Y ``=`` ``example_data``$``Y``,`\
`  g_learners ``=`` ``sl_lrnr``,`\
`  h_learners ``=`` ``sl_lrnr``,`\
`  b_learners ``=`` ``sl_lrnr``,`\
`  q_learners ``=`` ``sl_lrnr``,`\
`  r_learners ``=`` ``sl_lrnr``,`\
`  effect ``=`` ``"indirect"``,`\
`  estimator ``=`` ``"tmle"`\
`)`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``tmle_ie``)`

    ## # A tibble: 1 × 7
    ##   lwr_ci param_est  upr_ci var_est eif_mean estimator param                  
    ##    <dbl>     <dbl>   <dbl>   <dbl>    <dbl> <chr>     <chr>                  
    ## 1 -0.215    -0.134 -0.0528 0.00172   0.0213 tmle      indirect_interventional

From the output of the summary method, we note that the TML estimate of
the interventional indirect effect
$`\hat{\theta}_{\text{tmle}}^{\text{IE}}`$ is -0.134, with 95%
confidence interval \[-0.215, -0.053\]. As before, the TML estimator
provides better finite-sample performance than the one-step estimator,
so it may be preferred in this example.

### References

Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, et al. 2018.
“Double/Debiased Machine Learning for Treatment and Structural
Parameters.” *The Econometrics Journal* 21 (1).
<https://doi.org/10.1111/ectj.12097>.

Dı́az, Iván, Nima S Hejazi, Kara E Rudolph, and Mark J van der Laan.
2020. “Non-Parametric Efficient Causal Mediation with Intermediate
Confounders.” *Biometrika* 108 (3): 627–41.
<https://doi.org/10.1093/biomet/asaa085>.

Klaassen, Chris AJ. 1987. “Consistent Estimation of the Influence
Function of Locally Asymptotically Linear Estimators.” *Annals of
Statistics*, 1548–62. <https://doi.org/10.1214/aos/1176350609>.

Rudolph, Kara E, Oleg Sofrygin, Wenjing Zheng, and Mark J van der Laan.
2017. “Robust and Flexible Estimation of Stochastic Mediation Effects: A
Proposed Method and Example in a Randomized Trial Setting.”
*Epidemiologic Methods* 7 (1). <https://doi.org/10.1515/em-2017-0007>.

van der Laan, Mark J, Eric C Polley, and Alan E Hubbard. 2007. “Super
Learner.” *Statistical Applications in Genetics and Molecular Biology* 6
(1).

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

van der Laan, Mark J, and Sherri Rose. 2018. *Targeted Learning in Data
Science: Causal Inference for Complex Longitudinal Studies*. Springer
Science & Business Media.

VanderWeele, Tyler J, Stijn Vansteelandt, and James M Robins. 2014.
“Effect Decomposition in the Presence of an Exposure-Induced
Mediator-Outcome Confounder.” *Epidemiology* 25 (2): 300.
<https://doi.org/10.1097/ede.0000000000000034>.

Zheng, Wenjing, and Mark J van der Laan. 2011. “Cross-Validated Targeted
Minimum-Loss-Based Estimation.” In *Targeted Learning: Causal Inference
for Observational and Experimental Data*. Springer.
<https://doi.org/10.1007/978-1-4419-9782-1_27>.
