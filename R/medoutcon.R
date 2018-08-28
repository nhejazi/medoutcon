#' Nonparametric estimation of decomposition term for causal mediation analysis
#' with stochastic interventions under mediator-outcome confounding
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar corresponding to a
#'  set of baseline covariates.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param L A \code{numeric} vector corresponding to a mediator-outcome
#'  confounder affected by treatment (on the causal pathway between intervention
#'  A, mediator Z, and outcome Y, but unaffected itself by the mediator Z).
#' @param Z A \code{numeric} vector, \code{matrix}, \code{data.frame}, or
#'  similar corresponding to a set of mediators (on the causal pathway between
#'  the intervention A and the outcome Y).
#' @param Y A \code{numeric} vector corresponding to an outcome variable.
#' @param contrast ...
#' @param g_lrnrs A \code{Stack} object, or other learner class (inheriting from
#'  \code{Lrnr_base}), containing a single or set of instantiated learners from
#'  the \code{sl3} package, to be used in fitting a model for the propensity
#'  score, i.e., g = P(A | W).
#' @param e_lrnrs A \code{Stack} object, or other learner class (inheriting from
#'  \code{Lrnr_base}), containing a single or set of instantiated learners from
#'  the \code{sl3} package, to be used in fitting a cleverly parameterized
#'  propensity score that includes the mediators, i.e., e = P(A | Z, W).
#' @param q_lrnrs A \code{Stack} object, or other learner class (inheriting from
#'  \code{Lrnr_base}), containing a single or set of instantiated learners from
#'  the \code{sl3} package, to be used in fitting a regression involving the
#'  mediator-outcome confounder, i.e., q(L | A', W).
#' @param r_lrnrs A \code{Stack} object, or other learner class (inheriting from
#'  \code{Lrnr_base}), containing a single or set of instantiated learners from
#'  the \code{sl3} package, to be used in fitting a regression involving the
#'  mediator-outcome confounder, i.e., r(L | A', M, W).
#' @param u_lrnrs A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a reduced regression
#'  useful for computing the efficient one-step estimator, i.e., u(L, A, W) =
#'  E[m(A, L, Z, W) * (q(L|A,W) / r(L|A,Z,W)) * (e(a'|Z,W) / e(A|Z,W)) |
#'  L = l, A = a, W = w].
#' @param v_lrnrs A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a reduced regression
#'  useful for computing the efficient one-step estimator, i.e., v(A,W) =
#'  E[\int_z m(a', l, Z, W) * q(l|A',W) d\nu(z) | A = a, W = w)].
#' @param estimator The desired estimator of the natural direct effect to be
#'  computed. Currently, choices are limited to a substitution estimator, a
#'  re-weighted estimator, and an efficient one-step estimator. The interested
#'  user should consider consulting Díaz & Hejazi (2019+) for a comparative
#'  investigation of each of these estimators.
#' @param estimator_args A \code{list} of extra arguments to be passed (via
#'  \code{...}) to the function call for the specified estimator. The default
#'  is so chosen as to allow the number of folds used in computing the AIPW
#'  estimator to be easily tweaked. Refer to the documentation for functions
#'  \code{\link{est_onestep}}, \code{\link{est_ipw}}, and
#'  \code{\link{est_substitution}} for details on what other arguments may be
#'  specified through this mechanism.
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table as.data.table setnames
#' @importFrom origami make_folds cross_validate
#' @importFrom sl3 Lrnr_glm_fast
#' @importFrom stats binomial
#'
#' @export
#
medoutcon <- function(W,
                      A,
                      L,
                      Z,
                      Y,
                      contrast = c(0, 1),
                      g_lrnrs =
                        sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                      e_lrnrs =
                        sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                      q_lrnrs =
                        sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                      r_lrnrs =
                        sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                      u_lrnrs = sl3::Lrnr_glm_fast$new(),
                      v_lrnrs = sl3::Lrnr_glm_fast$new(),
                      estimator = c(
                        "onestep",
                        "tmle",
                        "sub",
                        "ipw"
                      ),
                      estimator_args = list(cv_folds = 5)) {
  # set defaults
  estimator <- match.arg(estimator)
  estimator_args <- unlist(estimator_args, recursive = FALSE)

  # NOTE: mediator-outcome confounding case
  # construct input data structure
  data <- data.table::as.data.table(cbind(Y, Z, L, A, W))
  w_names <- paste("W", seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  z_names <- paste("Z", seq_len(dim(data.table::as.data.table(Z))[2]),
    sep = "_"
  )
  data.table::setnames(data, c("Y", z_names, "L", "A", w_names))

  browser()

  if (estimator == "sub") {
    # SUBSTITUTION ESTIMATOR
    stop("The substitution estimator is currently under development.")
  } else if (estimator == "ipw") {
    # INVERSE PROBABILITY RE-WEIGHTED ESTIMATOR
    stop("The IPW estimator is currently under development.")
  } else if (estimator == "onestep") {
    # EFFICIENT ONE-STEP ESTIMATOR
    onestep_est_args <- list(
      data = data, contrast = contrast,
      g_lrnrs = g_lrnrs, e_lrnrs = e_lrnrs,
      q_lrnrs = q_lrnrs, r_lrnrs = r_lrnrs,
      u_lrnrs = u_lrnrs, v_lrnrs = v_lrnrs,
      w_names = w_names, z_names = z_names,
      estimator_args
    )
    est_out <- do.call(est_onestep_moc, onestep_est_args)
  } else if (estimator == "tmle") {
    # TARGETED MAXIMUM LIKELIHOOD ESTIMATOR (just call tmle3::fit_tmle?)
    stop("The TML estimator is currently under development.")
  }

  # lazily create output as classed list
  class(est_out) <- "medoutcon"
  return(est_out)
}
