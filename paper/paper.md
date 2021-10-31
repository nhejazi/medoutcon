---
title: "`medoutcon`: Efficient causal mediation analysis with machine learning in `R`"
tags:
  - causal inference
  - machine learning
  - semiparametric estimation
  - mediation analysis
  - natural direct effect
  - interventional direct effect
  - R
authors:
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1
  - name: Kara E. Rudolph
    orcid: 0000-0002-9417-7960
    affiliation: 2
  - name: Iván Díaz
    orcid: 0000-0001-9056-2047
    affiliation: 1
affiliations:
  - name: Division of Biostatistics, Department of Population Health Sciences, Weill Cornell Medicine
    index: 1
  - name: Department of Epidemiology, Mailman School of Public Health, Columbia University
    index: 2
date: 04 November 2021
bibliography: ../inst/REFERENCES.bib
---

# Summary

Efficient estimators of interventional (in)direct effects in the presence of
mediator-outcome confounding affected by exposure. The effects estimated allow
for the impact of the exposure on the outcome through a direct path to be
disentangled from that through mediators, even in the presence of intermediate
confounders that complicate such a relationship. Currently supported are
non-parametric efficient one-step and targeted minimum loss estimators based on
the formulation of Díaz, Hejazi, Rudolph, and van der Laan (2020)
<doi:10.1093/biomet/asaa085>. Support for efficient estimation of the natural
(in)direct effects is also provided, appropriate for settings in which
intermediate confounders are absent.

Statistical causal inference has traditionally focused on effects defined by
inflexible static interventions, applicable only to binary or categorical
exposures. The evaluation of such interventions is often plagued by many
problems, both theoretical (e.g., non-identification) and practical (e.g.,
positivity violations); however, stochastic interventions provide a promising
solution to these fundamental issues [@diaz2018stochastic]. The `txshift` `R`
package provides researchers in (bio)statistics, epidemiology, health policy,
economics, and related disciplines with access to state-of-the-art statistical
methodology for evaluating the causal effects of stochastic shift interventions
on _continuous-valued_ exposures. `txshift` estimates the causal effects of
modified treatment policies (or "feasible interventions"), which take into
account the natural value of an exposure in assigning an intervention level. To
accommodate use in study designs incorporating outcome-dependent two-phase
sampling (e.g., case-control), the package provides two types of modern
corrections, both rooted in semiparametric theory, for constructing unbiased and
efficient estimates, despite the significant limitations induced by such
designs. Thus, `txshift` makes possible the estimation of the causal effects of
stochastic interventions in experimental and observational study settings
subject to real-world design limitations that commonly arise in modern
scientific practice.

# Statement of Need

Researchers seeking to build upon or apply cutting-edge statistical approaches
for causal inference often face significant obstacles: such methods are usually
not accompanied by robust, well-tested, and well-documented software packages.
Yet coding such methods from scratch is often impractical for the applied
researcher, as understanding the theoretical underpinnings of these methods
requires advanced training, severely complicating the assessment and testing of
bespoke causal inference software. What's more, even when such software tools
exist, they are usually minimal implementations, providing support only for
deploying the statistical method in problem settings untouched by the
complexities of real-world data. The `txshift` `R` package solves this problem
by providing an open source tool for evaluating the causal effects of flexible,
stochastic interventions, applicable to categorical or continuous-valued
exposures, while providing corrections for appropriately handling data generated
by commonly used but complex two-phase sampling designs.

# Background

Causal inference has traditionally focused on the effects of static
interventions, under which the magnitude of the exposure is set to a fixed,
prespecified value for each unit. The evaluation of such interventions faces
a host of issues, among them non-identification, violations of the assumption of
positivity, and inefficiency. Stochastic interventions provide a promising
solution to these fundamental issues by allowing for the target parameter to be
defined as the mean counterfactual outcome under a hypothetically shifted
version of the observed exposure distribution [@diaz2012population].
Modified treatment policies, a particular class of such interventions, may be
interpreted as shifting the natural exposure level at the level of a given
observational unit [@haneuse2013estimation; @diaz2018stochastic].

Despite the promise of such advances in causal inference, real data analyses are
often further complicated by economic constraints, such as when the primary
variable of interest is far more expensive to collect than auxiliary covariates.
Two-phase sampling is often used to bypass these limitations -- unfortunately,
these sampling schemes produce side effects that require further adjustment when
formal statistical inference is the principal goal of a study. Among the rich
literature on two-phase designs, @rose2011targeted2sd stand out for providing
a study of nonparametric efficiency theory under a broad class of two-phase
designs. Their work provides guidance on constructing efficient estimators of
causal effects under general two-phase sampling designs.

# `txshift`'s Scope

Building on these prior works, @hejazi2020efficient outlined a novel approach
for use in such settings: augmented targeted minimum loss (TML) and one-step
estimators for the causal effects of stochastic interventions, with guarantees
of consistency, efficiency, and multiple robustness despite the presence of
two-phase sampling. These authors further outlined a technique that summarizes
the effect of shifting an exposure variable on the outcome of interest via
a nonparametric working marginal structural model, analogous to a dose-response
analysis. The `txshift` software package, for the `R` language and environment
for statistical computing [@R], implements this methodology.

`txshift` is designed to facilitate the construction of TML and one-step
estimators of the causal effects of modified treatment policies that shift the
observed exposure value up (or down) by an arbitrary scalar $\delta$, which may
possibly take into account the natural value of the exposure (and, in future
versions, the covariates). The `R` package includes tools for deploying these
efficient estimators under outcome-dependent two-phase sampling designs, with
two types of corrections: (1) a reweighting procedure that introduces inverse
probability of censoring weights directly into relevant loss functions, as
discussed in @rose2011targeted2sd; as well as (2) an augmented efficient
influence function estimating equation, studied more thoroughly by
@hejazi2020efficient. `txshift` integrates with the [`sl3`
package](https://github.com/tlverse/sl3) [@coyle2020sl3] to allow for ensemble
machine learning to be leveraged in the estimation of nuisance parameters.
What's more, the `txshift` package draws on both the `hal9001`
[@coyle2020hal9001; @hejazi2020hal9001] and `haldensify` [@hejazi2020haldensify]
`R` packages to allow each of the efficient estimators to be constructed in
a manner consistent with the methodological and theoretical advances of
@hejazi2020efficient, which require fast convergence rates of nuisance
parameters to their true counterparts for efficiency of the resultant estimator.

# Availability

The `txshift` package has been made publicly available both [via
GitHub](https://github.com/nhejazi/txshift) and the [Comprehensive `R` Archive
Network](https://CRAN.R-project.org/package=txshift). Use of the `txshift`
package has been extensively documented in the package's `README`, two
vignettes, and its [`pkgdown` documentation
website](https://code.nimahejazi.org/txshift).

# Acknowledgments

Nima Hejazi's contributions to this work were supported in part by a grant from
the National Institutes of Health: [T32
LM012417-02](https://projectreporter.nih.gov/project_info_description.cfm?aid=9248418&icde=37849831&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=).

# References

