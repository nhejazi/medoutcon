
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`medoutcon`

[![Travis-CI Build
Status](https://travis-ci.org/nhejazi/medoutcon.svg?branch=master)](https://travis-ci.org/nhejazi/medoutcon)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/medoutcon?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/medoutcon)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/medoutcon/master.svg)](https://codecov.io/github/nhejazi/medoutcon?branch=master)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Efficient Causal Mediation Analysis with Intermediate Confounders

**Authors:** [Nima Hejazi](https://nimahejazi.org), [Iván
Díaz](https://idiaz.xyz), and [Kara
Rudolph](http://biostat.jhsph.edu/~krudolph/)

-----

## What’s `medoutcon`?

The `medoutcon` R package provides facilities for efficient estimation
of stochastic (in)direct effects that measure the impact of a treatment
variable \(A\) on an outcome variable \(Y\), through a direct path
(through \(A\) only) and an indirect path (through a set of mediators
\(M\) only), in the presence of an intermediate
<b>med</b>iator-<b>out</b>come <b>con</b>founder \(Z\), itself affected
by the treatment \(A\). While the proposed approach is similar to those
appearing in VanderWeele, Vansteelandt, and Robins (2014), Rudolph et
al. (2017), and Zheng and van der Laan (2017), `medoutcon` is designed
as a software implementation to accompany the methodology proposed in
Díaz et al. (2019). Both an efficient one-step bias-corrected estimator
with cross-fitting (Pfanzagl and Wefelmeyer 1985; Zheng and van der Laan
2011; Chernozhukov et al. 2018) and a one-step cross-validated targeted
minimum loss (TML) estimator (van der Laan and Rose 2011; Zheng and van
der Laan 2011) are made available. `medoutcon` integrates with the
[`sl3` R package](https://github.com/tlverse/sl3) (Coyle et al. 2019) to
leverage statistical machine learning in the estimation procedure.

-----

## Installation

Install the most recent *stable release* from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("nhejazi/medoutcon")
```

-----

## Example

To illustrate how `medoutcon` may be used to estimate stochastic
(in)direct effects of the treatment (`A`) on the outcome (`Y`) in the
presence of mediator(s) (`M`) and a binary mediator-outcome confounder
(`Z`), consider the following simple example:

``` r
library(data.table)
library(tidyverse)
library(medoutcon)

# produces a simple data set based on ca causal model with mediation
make_example_data <- function(n_obs = 1000) {
  # baseline covariate -- simple, binary
  W <- replicate(2, rbinom(n_obs, 1, prob = 0.50))
  W <- as.data.table(W)
  setnames(W, c("w_1", "w_2"))

  # create treatment based on baseline W
  A <- as.numeric(rbinom(n_obs, 1, prob = (rowSums(W) / 3) + 0.1))

  # single mediator-outcome confounder
  z_prob <- 1 - plogis((A^2 + rowMeans(W)) / (A + rowSums(W^3) + 0.5))
  Z <- rbinom(n_obs, 1, prob = z_prob)

  # matrix of mediators
  m1_prob <- plogis((A^2 - Z) / (A + rowMeans(W) + 0.5))
  m2_prob <- 1 - plogis((A + rowSums(W)) / (Z + 0.5))
  M1 <- rbinom(n_obs, 1, prob = m1_prob)
  M2 <- rbinom(n_obs, 1, prob = m2_prob)
  M <- as.data.table(list(m_1 = M1, m_2 = M2))

  # create outcome as a linear function + white noise
  Y <- rowMeans(M) - Z + A - 0.1 * rowSums(W) +
    rnorm(n_obs, mean = 0, sd = 0.25)

  # full data structure
  data <- as.data.table(cbind(Y, M, Z, A, W))
  return(data)
}

# set seed and simulate example data
set.seed(75681)
example_data <- make_example_data(10000)
w_names <- str_subset(colnames(example_data), "w")
m_names <- str_subset(colnames(example_data), "m")

# compute one-step estimate
os_medoutcon <- medoutcon(W = example_data[, ..w_names],
                          A = example_data$A,
                          Z = example_data$Z,
                          M = example_data[, ..m_names],
                          Y = example_data$Y,
                          contrast = c(0, 1),
                          estimator = "onestep",
                          estimator_args = list(cv_folds = 3))
summary(os_medoutcon)
#>        lwr_ci     param_est        upr_ci     param_var      eif_mean 
#>       -0.2191       -0.1989       -0.1786         1e-04    1.4501e-17 
#>     estimator         param 
#>       onestep contrast_spec

# compute targeted minimum loss estimate
tmle_medoutcon <- medoutcon(W = example_data[, ..w_names],
                            A = example_data$A,
                            Z = example_data$Z,
                            M = example_data[, ..m_names],
                            Y = example_data$Y,
                            contrast = c(0, 1),
                            estimator = "tmle",
                            estimator_args = list(cv_folds = 3))
summary(tmle_medoutcon)
#>        lwr_ci     param_est        upr_ci     param_var      eif_mean 
#>       -0.3097       -0.2695       -0.2292         4e-04    3.7347e-18 
#>     estimator         param 
#>          tmle contrast_spec
```

For details on how to use data adaptive regression (machine learning)
techniques in the estimation of nuisance parameters, consider consulting
the vignette that accompanies the package.

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/medoutcon/issues).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/nhejazi/medoutcon/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `medoutcon` R package, please cite the following:

``` 
    @article{diaz2019nonparametric,
      title={Non-parametric efficient causal mediation with intermediate
        confounders},
      author={D{\'\i}az, Iv{\'a}n and Hejazi, Nima S and Rudolph, Kara E
        and {van der Laan}, Mark J},
      year={2019},
      url = {https://arxiv.org/abs/1912.09936},
      doi = {},
      journal={},
      volume={},
      number={},
      pages={},
      publisher={}
    }

    @manual{hejazi2020medoutcon,
      author={Hejazi, Nima S and D{\'\i}az, Iv{\'a}n and Rudolph, Kara E},
      title = {{medoutcon}: Efficient causal mediation analysis under
        intermediate confounding},
      year  = {2020},
      url = {https://github.com/nhejazi/medoutcon},
      note = {R package version 0.0.7}
    }
```

-----

## License

© 2019-2020 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2019-2020 Nima S. Hejazi
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

-----

## References

<div id="refs" class="references">

<div id="ref-chernozhukov2018double">

Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo,
Christian Hansen, Whitney Newey, and James Robins. 2018.
“Double/Debiased Machine Learning for Treatment and Structural
Parameters.” *The Econometrics Journal* 21 (1).
<https://doi.org/10.1111/ectj.12097>.

</div>

<div id="ref-coyle2019sl3">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, and Oleg Sofrygin. 2019.
“sl3: Modern Pipelines for Machine Learning and Super Learning.”
<https://github.com/tlverse/sl3>.
<https://doi.org/10.5281/zenodo.3558317>.

</div>

<div id="ref-diaz2019nonparametric">

Díaz, Iván, Nima S Hejazi, Kara E Rudolph, and Mark J van der Laan.
2019. “Non-Parametric Efficient Causal Mediation with Intermediate
Confounders.” <https://arxiv.org/abs/1912.09936>.

</div>

<div id="ref-pfanzagl1985contributions">

Pfanzagl, J, and W Wefelmeyer. 1985. “Contributions to a General
Asymptotic Statistical Theory.” *Statistics & Risk Modeling* 3 (3-4):
379–88.

</div>

<div id="ref-rudolph2017robust">

Rudolph, Kara E, Oleg Sofrygin, Wenjing Zheng, and Mark J van der Laan.
2017. “Robust and Flexible Estimation of Stochastic Mediation Effects: A
Proposed Method and Example in a Randomized Trial Setting.”
*Epidemiologic Methods* 7 (1). De Gruyter.

</div>

<div id="ref-vdl2011targeted">

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

</div>

<div id="ref-vanderweele2014effect">

VanderWeele, Tyler J, Stijn Vansteelandt, and James M Robins. 2014.
“Effect Decomposition in the Presence of an Exposure-Induced
Mediator-Outcome Confounder.” *Epidemiology (Cambridge, Mass.)* 25 (2).
NIH Public Access: 300.

</div>

<div id="ref-zheng2011cross">

Zheng, Wenjing, and Mark J van der Laan. 2011. “Cross-Validated Targeted
Minimum-Loss-Based Estimation.” In *Targeted Learning*, 459–74.
Springer.

</div>

<div id="ref-zheng2017longitudinal">

———. 2017. “Longitudinal Mediation Analysis with Time-Varying Mediators
and Exposures, with Application to Survival Outcomes.” *Journal of
Causal Inference* 5 (2). De Gruyter.

</div>

</div>
