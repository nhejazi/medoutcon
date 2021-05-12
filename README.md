
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`medoutcon`

[![Travis-CI Build
Status](https://travis-ci.com/nhejazi/medoutcon.svg?branch=master)](https://travis-ci.com/nhejazi/medoutcon)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/medoutcon?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/medoutcon)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/medoutcon/master.svg)](https://codecov.io/github/nhejazi/medoutcon?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Efficient Causal Mediation Analysis with Intermediate Confounders

**Authors:** [Nima Hejazi](https://nimahejazi.org), [Iván
Díaz](https://idiaz.xyz), and [Kara
Rudolph](https://kararudolph.github.io/)

-----

## What’s `medoutcon`?

The `medoutcon` R package provides facilities for efficient estimation
of path-specific (in)direct effects that measure the impact of a
treatment variable \(A\) on an outcome variable \(Y\), through a direct
path (through \(A\) only) and an indirect path (through a set of
mediators \(M\) only). In the presence of an intermediate
<b>med</b>iator-<b>out</b>come <b>con</b>founder \(Z\), itself affected
by the treatment \(A\), these correspond to the *interventional*
(in)direct effects described by Dı́az et al. (2020), though similar (yet
less general) effect definitions and/or estimation strategies have
appeared in VanderWeele, Vansteelandt, and Robins (2014), Rudolph et al.
(2017), Zheng and van der Laan (2017), and Benkeser (2020). When no
intermediate confounders are present, these effect definitions simplify
to the well-studied *natural* (in)direct effects, and our estimators are
analogs of those formulated by Zheng and van der Laan (2012). Both an
efficient one-step bias-corrected estimator with cross-fitting (Pfanzagl
and Wefelmeyer 1985; Zheng and van der Laan 2011; Chernozhukov et al.
2018) and a cross-validated targeted minimum loss estimator (TMLE) (van
der Laan and Rose 2011; Zheng and van der Laan 2011) are made available.
`medoutcon` integrates with the [`sl3` R
package](https://github.com/tlverse/sl3) (Coyle et al. 2020) to leverage
statistical machine learning in the estimation procedure.

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
example_data <- make_example_data()
w_names <- str_subset(colnames(example_data), "W")
m_names <- str_subset(colnames(example_data), "M")

# quick look at the data
head(example_data)
#>    W_1 W_2 W_3 A Z M Y
#> 1:   1   0   1 0 0 0 1
#> 2:   0   1   0 0 0 1 0
#> 3:   1   1   1 1 0 1 1
#> 4:   0   1   1 0 0 1 0
#> 5:   0   0   0 0 0 1 1
#> 6:   1   0   1 1 0 1 0

# compute one-step estimate of the interventional direct effect
os_de <- medoutcon(W = example_data[, ..w_names],
                   A = example_data$A,
                   Z = example_data$Z,
                   M = example_data[, ..m_names],
                   Y = example_data$Y,
                   effect = "direct",
                   estimator = "onestep")
os_de
#> $theta
#> [1] -0.07567883
#> 
#> $var
#> [1] 0.003171717
#> 
#> $type
#> [1] "onestep"
#> 
#> $param
#> [1] "direct_effect"

# compute targeted minimum loss estimate of the interventional direct effect
tmle_de <- medoutcon(W = example_data[, ..w_names],
                     A = example_data$A,
                     Z = example_data$Z,
                     M = example_data[, ..m_names],
                     Y = example_data$Y,
                     effect = "direct",
                     estimator = "tmle")
tmle_de
#> $theta
#> [1] -0.0775445
#> 
#> $var
#> [1] 0.003476221
#> 
#> $type
#> [1] "tmle"
#> 
#> $param
#> [1] "direct_effect"
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
    @article{diaz2020nonparametric,
      title={Non-parametric efficient causal mediation with intermediate
        confounders},
      author={D{\'\i}az, Iv{\'a}n and Hejazi, Nima S and Rudolph, Kara E
        and {van der Laan}, Mark J},
      year={2020},
      url = {https://arxiv.org/abs/1912.09936},
      doi = {10.1093/biomet/asaa085},
      journal={Biometrika},
      volume={},
      number={},
      pages={},
      publisher={Oxford University Press}
    }

    @software{hejazi2021medoutcon,
      author={Hejazi, Nima S and D{\'\i}az, Iv{\'a}n and Rudolph, Kara E},
      title = {{medoutcon}: Efficient causal mediation analysis under
        intermediate confounding},
      year  = {2021},
      url = {https://github.com/nhejazi/medoutcon},
      note = {R package version 0.1.5}
    }
```

-----

## License

© 2020-2021 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2020-2021 Nima S. Hejazi
    
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

<div id="ref-benkeser2020nonparametric">

Benkeser, David. 2020. “Nonparametric Inference for Interventional
Effects with Multiple Mediators.” *arXiv Preprint arXiv:2001.06027*.

</div>

<div id="ref-chernozhukov2018double">

Chernozhukov, Victor, Denis Chetverikov, Mert Demirer, Esther Duflo,
Christian Hansen, Whitney Newey, and James Robins. 2018.
“Double/Debiased Machine Learning for Treatment and Structural
Parameters.” *The Econometrics Journal* 21 (1).
<https://doi.org/10.1111/ectj.12097>.

</div>

<div id="ref-coyle2020sl3">

Coyle, Jeremy R, Nima S Hejazi, Ivana Malenica, and Oleg Sofrygin. 2020.
“`sl3`: Modern Pipelines for Machine Learning and Super Learning.”
<https://github.com/tlverse/sl3>.
<https://doi.org/10.5281/zenodo.1342293>.

</div>

<div id="ref-diaz2020nonparametric">

Dı́az, Iván, Nima S Hejazi, Kara E Rudolph, and Mark J van der Laan.
2020. “Non-Parametric Efficient Causal Mediation with Intermediate
Confounders.” *Biometrika*. <https://doi.org/10.1093/biomet/asaa085>.

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
*Epidemiologic Methods* 7 (1).

</div>

<div id="ref-vdl2011targeted">

van der Laan, Mark J, and Sherri Rose. 2011. *Targeted Learning: Causal
Inference for Observational and Experimental Data*. Springer Science &
Business Media.

</div>

<div id="ref-vanderweele2014effect">

VanderWeele, Tyler J, Stijn Vansteelandt, and James M Robins. 2014.
“Effect Decomposition in the Presence of an Exposure-Induced
Mediator-Outcome Confounder.” *Epidemiology (Cambridge, Mass.)* 25 (2):
300.

</div>

<div id="ref-zheng2011cross">

Zheng, Wenjing, and Mark J van der Laan. 2011. “Cross-Validated Targeted
Minimum-Loss-Based Estimation.” In *Targeted Learning*, 459–74.
Springer.

</div>

<div id="ref-zheng2012targeted">

———. 2012. “Targeted Maximum Likelihood Estimation of Natural Direct
Effects.” *International Journal of Biostatistics* 8 (1).
<https://doi.org/10.2202/1557-4679.1361>.

</div>

<div id="ref-zheng2017longitudinal">

———. 2017. “Longitudinal Mediation Analysis with Time-Varying Mediators
and Exposures, with Application to Survival Outcomes.” *Journal of
Causal Inference* 5 (2).

</div>

</div>
