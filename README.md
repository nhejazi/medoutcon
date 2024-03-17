
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`medoutcon`

<!-- badges: start -->

[![R-CMD-check](https://github.com/nhejazi/medoutcon/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/nhejazi/medoutcon/actions/workflows/R-CMD-check.yml)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/medoutcon/master.svg)](https://codecov.io/github/nhejazi/medoutcon?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5809519.svg)](https://doi.org/10.5281/zenodo.5809519)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03979/status.svg)](https://doi.org/10.21105/joss.03979)
<!-- badges: end -->

> Efficient Causal Mediation Analysis for the Natural and Interventional
> Effects

**Authors:** [Nima Hejazi](https://nimahejazi.org), [Iván
Díaz](https://idiaz.xyz), and [Kara
Rudolph](https://kararudolph.github.io/)

------------------------------------------------------------------------

## What’s `medoutcon`?

The `medoutcon` R package provides facilities for efficient estimation
of path-specific (in)direct effects that measure the impact of a
treatment variable $A$ on an outcome variable $Y$, through a direct path
(through $A$ only) and an indirect path (through a set of mediators $M$
only). In the presence of an intermediate <b>med</b>iator-<b>out</b>come
<b>con</b>founder $Z$, itself affected by the treatment $A$, these
correspond to the *interventional* (in)direct effects described by Dı́az
et al. (2020), though similar (yet less general) effect definitions
and/or estimation strategies have appeared i`n @`vanderweele2014effect,
Rudolph et al. (2017), Zheng and van der Laan (2017), and Benkeser and
Ran (2021). When no intermediate confounders are present, these effect
definitions simplify to the well-studied *natural* (in)direct effects,
and our estimators are analogs of those formulated by Zheng and van der
Laan (2012). Both an efficient one-step bias-corrected estimator with
cross-fitting (Pfanzagl and Wefelmeyer 1985; Zheng and van der Laan
2011; Chernozhukov et al. 2018) and a cross-validated targeted minimum
loss estimator (TMLE) (van der Laan and Rose 2011; Zheng and van der
Laan 2011) are made available. `medoutcon` integrates with the [`sl3` R
package](https://github.com/tlverse/sl3) (Coyle et al. 2021) to leverage
statistical machine learning in the estimation procedure.

------------------------------------------------------------------------

## Installation

Install the most recent *stable release* from GitHub via
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("nhejazi/medoutcon")
```

------------------------------------------------------------------------

## Example

To illustrate how `medoutcon` may be used to estimate stochastic
interventional (in)direct effects of the exposure (`A`) on the outcome
(`Y`) in the presence of mediator(s) (`M`) and a mediator-outcome
confounder (`Z`), consider the following example:

``` r
library(data.table)
library(tidyverse)
library(medoutcon)
set.seed(1584)

# produces a simple data set based on ca causal model with mediation
make_example_data <- function(n_obs = 1000) {
  ## baseline covariates
  w_1 <- rbinom(n_obs, 1, prob = 0.6)
  w_2 <- rbinom(n_obs, 1, prob = 0.3)
  w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
  w <- cbind(w_1, w_2, w_3)
  w_names <- paste("W", seq_len(ncol(w)), sep = "_")

  ## exposure
  a <- as.numeric(rbinom(n_obs, 1, plogis(rowSums(w) - 2)))

  ## mediator-outcome confounder affected by treatment
  z <- rbinom(n_obs, 1, plogis(rowMeans(-log(2) + w - a) + 0.2))

  ## mediator -- could be multivariate
  m <- rbinom(n_obs, 1, plogis(rowSums(log(3) * w[, -3] + a - z)))
  m_names <- "M"

  ## outcome
  y <- rbinom(n_obs, 1, plogis(1 / (rowSums(w) - z + a + m)))

  ## construct output
  dat <- as.data.table(cbind(w = w, a = a, z = z, m = m, y = y))
  setnames(dat, c(w_names, "A", "Z", m_names, "Y"))
  return(dat)
}

# set seed and simulate example data
example_data <- make_example_data(n_obs = 5000L)
w_names <- str_subset(colnames(example_data), "W")
m_names <- str_subset(colnames(example_data), "M")

# quick look at the data
head(example_data)
#>    W_1 W_2 W_3 A Z M Y
#> 1:   1   1   0 0 1 0 1
#> 2:   0   0   0 1 0 1 0
#> 3:   1   0   1 0 0 1 0
#> 4:   0   0   0 0 0 1 0
#> 5:   0   0   0 0 0 0 1
#> 6:   1   0   0 1 0 1 0

# compute one-step estimate of the interventional direct effect
os_de <- medoutcon(
  W = example_data[, ..w_names],
  A = example_data$A,
  Z = example_data$Z,
  M = example_data[, ..m_names],
  Y = example_data$Y,
  effect = "direct",
  estimator = "onestep"
)
os_de

# compute targeted minimum loss estimate of the interventional direct effect
tmle_de <- medoutcon(
  W = example_data[, ..w_names],
  A = example_data$A,
  Z = example_data$Z,
  M = example_data[, ..m_names],
  Y = example_data$Y,
  effect = "direct",
  estimator = "tmle"
)
tmle_de
```

For details on how to use data adaptive regression (machine learning)
techniques in the estimation of nuisance parameters, consider consulting
the vignette that accompanies the package.

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/medoutcon/issues).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/nhejazi/medoutcon/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

------------------------------------------------------------------------

## Citation

After using the `medoutcon` R package, please cite the following:

        @article{diaz2020nonparametric,
          title={Non-parametric efficient causal mediation with intermediate
            confounders},
          author={D{\'\i}az, Iv{\'a}n and Hejazi, Nima S and Rudolph, Kara E
            and {van der Laan}, Mark J},
          year={2020},
          url = {https://arxiv.org/abs/1912.09936},
          doi = {10.1093/biomet/asaa085},
          journal={Biometrika},
          volume = {108},
          number = {3},
          pages = {627--641},
          publisher={Oxford University Press}
        }

        @article{hejazi2022medoutcon-joss,
          author = {Hejazi, Nima S and Rudolph, Kara E and D{\'\i}az,
            Iv{\'a}n},
          title = {{medoutcon}: Nonparametric efficient causal mediation
            analysis with machine learning in {R}},
          year = {2022},
          doi = {10.21105/joss.03979},
          url = {https://doi.org/10.21105/joss.03979},
          journal = {Journal of Open Source Software},
          publisher = {The Open Journal}
        }

        @software{hejazi2022medoutcon-rpkg,
          author={Hejazi, Nima S and D{\'\i}az, Iv{\'a}n and Rudolph, Kara E},
          title = {{medoutcon}: Efficient natural and interventional causal
            mediation analysis},
          year  = {2022},
          doi = {10.5281/zenodo.5809519},
          url = {https://github.com/nhejazi/medoutcon},
          note = {R package version 0.1.6}
        }

------------------------------------------------------------------------

## License

© 2020-2022 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License

    Copyright (c) 2020-2022 Nima S. Hejazi

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

------------------------------------------------------------------------

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-benkeser2020nonparametric" class="csl-entry">

Benkeser, David, and Jialu Ran. 2021. “Nonparametric Inference for
Interventional Effects with Multiple Mediators.” *Journal of Causal
Inference*. <https://doi.org/10.1515/jci-2020-0018>.

</div>

<div id="ref-chernozhukov2018double" class="csl-entry">

Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo,
Christian Hansen, Whitney Newey, and James Robins. 2018.
“Double/Debiased Machine Learning for Treatment and Structural
Parameters.” *The Econometrics Journal* 21 (1).
<https://doi.org/10.1111/ectj.12097>.

</div>

<div id="ref-coyle-gh-sl3" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, Rachael V Phillips, and
Oleg Sofrygin. 2021. “`sl3`: Modern Machine Learning Pipelines for Super
Learning.” <https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-diaz2020nonparametric" class="csl-entry">

Dı́az, Iván, Nima S Hejazi, Kara E Rudolph, and Mark J van der Laan.
2020. “Non-Parametric Efficient Causal Mediation with Intermediate
Confounders.” *Biometrika* 108 (3): 627–41.
<https://doi.org/10.1093/biomet/asaa085>.

</div>

<div id="ref-pfanzagl1985contributions" class="csl-entry">

Pfanzagl, J, and W Wefelmeyer. 1985. “Contributions to a General
Asymptotic Statistical Theory.” *Statistics & Risk Modeling* 3 (3-4):
379–88. <https://doi.org/10.1007/978-1-4612-5769-1>.

</div>

<div id="ref-rudolph2017robust" class="csl-entry">

Rudolph, Kara E, Oleg Sofrygin, Wenjing Zheng, and Mark J van der Laan.
2017. “Robust and Flexible Estimation of Stochastic Mediation Effects: A
Proposed Method and Example in a Randomized Trial Setting.”
*Epidemiologic Methods* 7 (1). <https://doi.org/10.1515/em-2017-0007>.

</div>

<div id="ref-vdl2011targeted" class="csl-entry">

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

</div>

<div id="ref-zheng2011cross" class="csl-entry">

Zheng, Wenjing, and Mark J van der Laan. 2011. “Cross-Validated Targeted
Minimum-Loss-Based Estimation.” In *Targeted Learning: Causal Inference
for Observational and Experimental Data*, 459–74. Springer.
<https://doi.org/10.1007/978-1-4419-9782-1_27>.

</div>

<div id="ref-zheng2012targeted" class="csl-entry">

———. 2012. “Targeted Maximum Likelihood Estimation of Natural Direct
Effects.” *International Journal of Biostatistics* 8 (1).
<https://doi.org/10.2202/1557-4679.1361>.

</div>

<div id="ref-zheng2017longitudinal" class="csl-entry">

———. 2017. “Longitudinal Mediation Analysis with Time-Varying Mediators
and Exposures, with Application to Survival Outcomes.” *Journal of
Causal Inference* 5 (2). <https://doi.org/10.1515/jci-2016-0006>.

</div>

</div>
