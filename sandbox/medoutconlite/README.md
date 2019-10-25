
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`medoutconlite`

> Efficient Causal Mediation Analysis with Mediator-Outcome Confounding

**Authors:** [Nima Hejazi](https://nimahejazi.org), [Iván
Díaz](https://idiaz.xyz), and [Kara
Rudolph](http://biostat.jhsph.edu/~krudolph/)

-----

## What’s `medoutconlite`?

The `medoutconlite` R package provides facilities for efficient
estimation of stochastic (in)direct effects that measure the impact of a
treatment variable \(A\) on an outcome variable \(Y\), through a direct
path (through \(A\) only) and an indirect path (through a set of
mediators \(M\) only), in the presence of an intermediate
<b>med</b>iator-<b>out</b>come <b>con</b>founder \(Z\), itself affected
by the treatment \(A\). It serves as a lite version of the [`medoutcon`
package](https://github.com/nhejazi/medoutcon) with fewer dependencies.

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
library(SuperLearner)
devtools::load_all(".")
#source('utils.r')
#source('slfunctions.r')
#source('estimators.r')
candidates <- c("SL.glm", "SL.caretRF", "SL.caretXGB", "SL.glmnet", "SL.earth")

n <- 500
data <- simdata(n)
effects <- mediation(data, weights = rep(1, n),
                     candidatesg = candidates,
                     candidatese = candidates,
                     candidatesm = candidates,
                     candidatesr = candidates,
                     candidatesq = candidates,
                     candidatesu = candidates,
                     candidatesv = candidates,
                     nfolds = 5,
                     family.outcome = binomial())




# produces a simple data set based on ca causal model with mediation
make_example_data <- function(n_obs = 1000) {
  # baseline covariate -- simple, binary
  W <- replicate(2, rbinom(n_obs, 1, prob = 0.50))
  W <- as.data.table(W)
  setnames(W, c("w_1", "w_2"))

  # create treatment based on baseline W
  A <- as.numeric(rbinom(n_obs, 1, prob = rowSums(W)/3 + 0.1))

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
```
