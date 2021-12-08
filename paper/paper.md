---
title: "`medoutcon`: Nonparametric efficient causal mediation analysis with machine learning in `R`"
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
applications, only the portion of the causal effect of an exposure on an outcome
through a particular pathway under study is of interest. The study of such
path-specific, or mediation, effects has a rich history, first undertaken
scientifically by @wright1921correlation and @wright1934method. Today, the study
of such effects has attracted a great deal of attention in statistics and causal
inference, inspired by applications in disciplines ranging from epidemiology and
vaccinology to psychology and economics. Examples include understanding the
biological mechanisms by which vaccines causally alter infection risk
[@hejazi2020efficient; @benkeser2021inference], assessing the effect of novel
pharmacological therapies on substance abuse disorder relapse
[@hejazi2021nonparametric; @rudolph2020explaining], and evaluating the effects
of housing vouchers on adolescent development [@rudolph2021helped]. The
`medoutcon` `R` package provides researchers in each of these disciplines, and
in others, with the tools necessary to implement statistically efficient
estimators of the _interventional_ direct and indirect effects
[@diaz2020nonparametric] (for brevity, henceforth, (in)direct effects),
a recently formulated set of causal effects robust to the presence of
confounding of the mediator-outcome relationship by the exposure. In cases where
such confounding is a nonissue, the interventional (in)direct effects
[@vanderweele2014effect] reduce to the well-studied _natural_ (in)direct effects
[@robins1992identifiability; @pearl2001direct], for which `medoutcon` provides
efficient estimators similar to those of @zheng2012targeted. By readily
incorporating the use of machine learning in the estimation of nuisance
parameters (through integration with the `sl3` `R` package [@coyle-gh-sl3] of
the `tlverse` ecosystem [@vdl2022targeted]), `medoutcon` incorporates
state-of-the-art non/semi-parametric estimation techniques, facilitating their
adoption in a vast array of settings.

# Statement of Need

While there is demonstrable interest in causal mediation analysis in a large
variety of disciplines, thoughtfully implementing data analysis strategies based
on recent developments in this area is challenging. Contributions in the causal
inference and statistics literature largely fall into two key areas. Broadly,
the study of identification outlines novel causal effect parameters with
properties desirable in real-world settings (e.g., the interventional effects,
which can be learned under mediator-outcome confounding) and untestable
assumptions under which a statistical functional corresponds to a causal
estimand of interest. A complementary line of study develops non/semi-parametric
efficiency theory for the statistical functionals outlined in the causal
identification literature, allowing for their robust estimation with modern
techniques from machine learning. Neither concerns itself with opening the door
to applying these estimators in real-world data analyses. Moreover, the
implementation of open source software for efficient estimators of causal
effects is complex -- for such a task, the data scientist must be
knowledgeable of causal inference, semiparametric statistical theory, machine
learning, and the intersection of these disciplines, and that is to forego
mention of research software engineering best practices, including, for example,
unit/regression testing and automated continuous integration. The `medoutcon`
`R` package is a free, open source implementation of non/semi-parametric
efficient estimators of the natural and interventional (in)direct effects,
providing data scientists in research and in industry with access to
state-of-the-art statistical methodology for causal mediation analysis. Its
estimators have been interrogated in simulation studies and applied in
real-world data analyses. To the best of our knowledge, no other `R` package
provides similarly convenient access to multiply robust, non/semi-parametric
efficient estimators of causal mediation effects with a flexible interface to
accommodate machine learning of nuisance parameters.

# Natural and Interventional Causal Mediation Effects

To evaluate the causal effects of an exposure on an outcome through mediating
pathways, let's consider a dataset of $n$ units, where the observed data on
a single unit is assumed to have been generated by a nonparametric structural
equation model (NPSEM) [@pearl2009causality]:
\begin{align*}
  W &= f_W(U_W); A = f_A(W, U_A); Z=f_Z(W, A, U_Z);\\
  M &= f_M(W, A, Z, U_M); Y = f_Y(W, A, Z, M, U_Y),
\end{align*}
where $W$ are baseline (pre-exposure) covariates, $A \in \{0,1\}$ is the
(binary) exposure of interest, $Z$ is an intermediate confounder of the
mediator-outcome relationship and is affected by exposure $A$, $M$ represents
mediating variables, and $Y$ is the outcome. This NPSEM admits an equivalent
representation as a directed acyclic graph (or DAG), in which each variable is
a node and dependencies are represented by directed paths between the nodes. The
natural (in)direct effects cannot generally be identified (i.e., learned from
the observed data) in the presence of intermediate confounding, so, for now, we
make the simplifying assumption that the intermediate variable $Z$ is absent. In
this simple case, the population average treatment effect (ATE) -- that is, the
total effect of $A$ on $Y$, comparing two exposure contrasts $\{a', a^{\star}\}$
-- may be decomposed into the natural direct effect (NDE) and the natural
indirect effect (NIE) as
\begin{equation*}
  \mathbb{E}[Y(a') - Y(a^{\star})] =
    \underbrace{\mathbb{E}[Y(a', M(a')) - Y(a',
      M(a^{\star}))]}_{\text{Indirect effect (through $M$)}} +
    \underbrace{\mathbb{E}[Y(a', M(a^{\star})) - Y(a^{\star},
      M(a^{\star}))]}_{\text{Direct effect (not through $M$)}},
\end{equation*}
where the _counterfactual_ variables $Y(\cdot)$ are _potential outcomes_
[@imbens2015causal; @hernan2021causal] -- that is, $Y(a')$ is the value that the
outcome would take when the exposure is set to level $a'$, possibly contrary to
fact. Similarly, $M(a^{\star})$ is the value that the mediators would take when
the exposure is set to level $a^{\star}$, as the result of an intervention, for
example. The NIE captures the effect of the exposure $A$ on $Y$ through the
mediating variables $M$ while the NDE captures the effect of $A$ on $Y$ through
all other pathways. @robins1992identifiability and @pearl2001direct
independently studied this decomposition within the potential outcomes and NPSEM
frameworks, respectively. In both cases, the NDE and NIE are derived from the
ATE by introducing a decomposition term that deterministically sets the values
of the exposure and mediators to differing values by the application of _static_
interventions. As regards estimation, @tchetgen2012semiparametric and
@zheng2012targeted outlined non/semi-parametric efficiency theory for developing
estimators of the NDE and NIE and proposed efficient estimators of these causal
quantities.

The presence of intermediate confounders $Z$ often cannot be ruled out in
real-world data analysis scenarios. Such post-exposure variables, which are
affected by $A$ and affect both $M$ and $Y$, complicate efforts to disentangle
the effect of $A$ on $Y$ through paths involving $M$ and other paths.
Recognizing the limitations of the natural effects in these settings,
@didelez2006direct, @petersen2006estimation, @vanderweele2014effect, and
@rudolph2017robust, among others, contributed to the development of the
interventional (in)direct effects. Unlike the decomposition strategy that
delineates the NDE and NIE, these effects require a more sophisticated approach
to identification, relying upon _stochastic_ interventions on the mediator(s),
which require random draws from the mediators post-intervention distribution
rather than the setting of fixed counterfactual values. Specifically, for the
two exposure contrasts $\{a', a^{\star}\}$, the effect of $A$ on $Y$ can be
defined as the difference in expected outcome in the hypothetical worlds in
which $(A,M) = (a', G_{a'})$ versus $(A,M) = (a^{\star}, G_{a^{\star}})$. Here,
$G_a$ denotes a random draw from the conditional distribution of $M_a$
conditional on $W$, as defined by a stochastic intervention. The direct and
indirect effects are defined as follows
\begin{equation*}
\mathbb{E}[Y(a', G_{a'}) - Y(a^{\star}, G_{a^{\star}})] =
  \underbrace{\mathbb{E}[Y(a', G_{a'}) - Y(a',
    G_{a^{\star}})]}_{\text{Indirect effect (through $M$)}} +
  \underbrace{\mathbb{E}[Y(a', G_{a^{\star}}) - Y(a^{\star},
      G_{a^{\star}})]}_{\text{Direct effect (not through $M$)}}.
\end{equation*}
Like the NDE, this interventional direct effect measures the effects through all
paths avoiding the mediating variables. Analogous to the NIE, the interventional
indirect effect measures the effect through paths involving the mediators. Note,
however, that natural and interventional mediation effects have different
interpretations. That is, the interventional indirect effect measures the effect
of fixing the exposure at $a'$ while setting the mediator to a random draw
$G_{a^{\star}}$ (i.e., under an intervention setting the exposure to
$a^{\star}$) versus a random draw $G_{a'}$ (i.e., after setting the exposure to
$a'$), given covariates $W$. Intuitively, the interventional effects remain
identifiable under intermediate confounding since the stochastic intervention on
the mediators breaks the relationship between $Z$ and $M$. Prior to the work of
@diaz2020nonparametric, and contemporaneous developments by
@benkeser2020nonparametric, non/semi-parametric efficiency theory for the
interventional (in)direct effects was unavailable. Recently, a novel family of
interventional effects, accommodating flexible stochastic interventions on the
exposure, have been formulated [@hejazi2021nonparametric].

# `medoutcon`'s Scope

Development of the `medoutcon` package began as a software accompaniment to the
theoretical developments of @diaz2020nonparametric -- where the investigations
of these authors outlined efficient estimators of the interventional (in)direct
effects, `medoutcon` implements these efficient estimators. Implemented in the
`R` language and environment for statistical computing [@R], `medoutcon` aims to
provide a simple application user interface (API) for convenience in a variety
of data analytic applications. Specifically, `medoutcon` -- via a single,
user-facing eponymous function `medoutcon()` -- provides access to both one-step
and targeted minimum loss (TML) estimators of these causal (in)direct effects.
State-of-the-art machine learning algorithms, including ensemble modeling
[@vdl2007super], may readily be used for the estimation of relevant nuisance
parameters, through a design that tightly couples `medoutcon` with the `sl3` `R`
package [@coyle-gh-sl3]. Cross-fitting is automatically incorporated, via the
`origami` `R` package [@coyle2018origami; @coyle-cran-origami], in computing the
efficient estimators, allowing for some common but restrictive regularity
conditions to be relaxed [@bickel1993efficient; @zheng2011cross;
@chernozhukov2017double].

Beyond implementing the interventional (in)direct effects, `medoutcon`
additionally allows for the natural (in)direct effects to be estimated when
intermediate confounders are omitted from the call to the `medoutcon()` function
(i.e., by setting `Z = NULL`). This feature is based on a correspondence between
the identifying statistical functionals of the natural and interventional
(in)direct effects in the absence of intermediate confounding. In this
simplified case, the efficient estimators of the interventional (in)direct
effects formulated by @diaz2020nonparametric are analogous to the efficient
estimators of the natural (in)direct effects formulated by @zheng2012targeted.
By supporting this case, `medoutcon` serves as a one-stop tool for estimating
these classical and popular causal mediation effects, allowing for practicing
data scientists and applied statisticians to deploy cutting-edge estimators of
the natural and interventional (in)direct effects through a unified API.

# Availability

The `medoutcon` package is publicly available [via
GitHub](https://github.com/nhejazi/medoutcon), with plans for submission to the
Comprehensive `R` Archive Network, pending the inclusion of its dependencies
(`sl3`, in particular) in that repository. Use of the `medoutcon` package has
been extensively documented in the package's `README`, a vignette, and its
[documentation website](https://code.nimahejazi.org/medoutcon). Ongoing
development of the package incorporates research and data science software
engineering best practices, including a suite of unit tests and automated
continuous integration checking.

# Acknowledgments

NSH's contributions to this work were supported in part by a grant from the
National Science Foundation (award number [DMS
2102840](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2102840)).

# References

