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
date: 07 November 2021
bibliography: ../inst/REFERENCES.bib
---

# Summary

Science is most often concerned with questions of _mechanism_. In myriad
applications, only the portion of the effect of an exposure variable on an
outcome variable through a particular pathway under study is of interest. The
study of such path-specific, or mediation, effects has a rich history, first
undertaken scientifically by wright1 and extended soon thereafter in wright2.
Today, the study of such effects has attracted a great deal of attention in
statistics and causal inference, often inspired by applications in disciplines
from epidemiology and vaccinology to psychology and economics. Examples include
understanding the biological mechanisms by which vaccines causally alter
infection risk [@hejazi2020efficient; benkeser], assessing the effect of novel
pharmacological therapies on substance abuse disorder relapse [@hejazi2021], and
evaluating the effects of housing vouchers on adolescent development [@rudolph].
The `medoutcon` `R` package provides researchers in each of these disciplines
with the tools necessary to implement efficient estimators of the interventional
(in)direct effects [diaz2020biometrika], a recently formulated set of causal
effects robust to the presence of confounding of the mediator-outcome
relationship by the exposure variable. In cases where such confounding is not
a problem, the interventional (in)direct effects [tyler, many others] reduce to
the well-studied natural (in)direct effects [robins, pearl], for which
`medoutcon` provides efficient estimators similar to those of @zheng. By readily
incorporating the use of machine learning in the estimation of nuisance
parameters (through integration with the `sl3` `R` package [coyle] of the
`tlverse` ecosystem [coyle]), `medoutcon` furnishes both researchers and
analysts with access to state-of-the-art semiparametric estimation techniques,
facilitating their use in a vast range of subject areas.

# Statement of Need

While there is demonstrable interest in causal mediation analysis in a large
variety of disciplines, thoughtfully implementing data analysis strategies based
on recent developments in this area is challenging. Developments in the causal
inference and statistics literature often fall into two key areas. Broadly, the
study of identification outlines novel causal effect parameters with properties
desirable in real-world settings (e.g., the interventional effects, which can be
learned under mediator-outcome confounding) and untestable assumptions under
which a statistical functionals corresponds to a given causal effect. Another
line of study develops semiparametric efficiency theory for the statistical
functionals matching these novel causal estimands, allowing their robust
estimation with tools from machine learning. While these complementary efforts
are necessary, neither is concerned with opening the door to applying these
estimators in real-world data analyses. What's more, the implementation of
open source software for efficient estimators of causal effects is no easy feat
-- for such a task, the data scientist must be knowledgeable of causal
inference, semiparametric theory, machine learning, and the intersection of
these disciplines (not to mention research software engineering best practices,
including, for example, unit/regression testing and continuous integration).
With these issues in mind, the `medoutcon` `R` package is a free, open source
implementation of non/semi-parametric efficient estimators of the natural and
interventional (in)direct effects, providing data scientists in research and in
industry with access to state-of-the-art statistical methodology for causal
mediation analysis.

# The Natural and Interventional Effects

Paragraph on the Natural effects

Paragraph on the Interventional Effects

# `medoutcon`'s Scope

Building on existing efforts in the literature on the interventional effects,
@diaz2020nonparametric recently developed non/semi-parametric efficiency theory

outlined a novel approach
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

The `medoutcon` package has been made publicly available [via
GitHub](https://github.com/nhejazi/medoutcon), with plans for submission to the
Comprehensive `R` Archive Network, pending the inclusion of its dependencies in
that repository. Use of the `medoutcon` package has been extensively documented
in the package's `README`, a vignette, and its [`pkgdown` documentation
website](https://code.nimahejazi.org/medoutcon).

# Acknowledgments

NH's contributions to this work were supported in part by a grant from the
National Science Foundation (award number DMS TODO).

# References

