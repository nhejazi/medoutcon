#' Defines a TML Estimator for Outcome under Joint Stochastic Intervention on
#' Treatment and Mediator
#'
#' @importFrom R6 R6Class
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'
#' @export
#
tmle3_Spec_medshift <- R6::R6Class(
  classname = "tmle3_Spec_medoutcon",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(shift_type = "exptilt", delta,
                              e_learners, phi_learners,
                              max_iter = 1e4, step_size = 1e-6,
                              ...) {
      options <- list(
        shift_type = shift_type,
        delta_shift = delta,
        e_learners = e_learners,
        phi_learners = phi_learners,
        max_iter = max_iter,
        step_size = step_size,
        ...
      )
      do.call(super$initialize, options)
    },
    make_tmle_task = function(data, node_list, ...) {
      # get variable types by guessing
      variable_types <- self$options$variable_types

      # build custom NPSEM including mediators with helper function
      npsem <- stochastic_mediation_npsem(node_list)

      # set up TMLE task based on NPSEM and return
      tmle_task <- tmle3_Task$new(data, npsem, variable_types)
      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # build likelihood using helper function and return
      likelihood <- stochastic_mediation_likelihood(tmle_task, learner_list)
      return(likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      # add derived likelihood factors to targeted likelihood object
      lf_e <- tmle3::define_lf(
        tmle3::LF_derived, "E", self$options$e_learners,
        targeted_likelihood, make_e_task
      )
      lf_phi <- tmle3::define_lf(
        tmle3::LF_derived, "phi", self$options$phi_learners,
        targeted_likelihood, make_phi_task
      )
      targeted_likelihood$add_factors(lf_e)
      targeted_likelihood$add_factors(lf_phi)

      # compute a tmle3 "by hand"
      tmle_params <- tmle3::define_param(Param_medoutcon, targeted_likelihood,
        shift_param = self$options$delta_shift
      )
      tmle_params <- list(tmle_params)
      return(tmle_params)
    },
    make_updater = function() {
      # default to ULFM approach
      updater <- tmle3_Update$new(
        one_dimensional = TRUE,
        constrain_step = TRUE,
        maxit = self$options$max_iter,
        delta_epsilon = self$options$step_size,
        cvtmle = TRUE
      )
    }
  ),
  active = list(),
  private = list()
)

################################################################################

#' Outcome under Joint Stochastic Intervention on Treatment and Mediator
#'
#' O = (W, A, Z, Y)
#' W = Covariates (possibly multivariate)
#' A = Treatment (binary or categorical)
#' Z = Mediators (binary or categorical; possibly multivariate)
#' Y = Outcome (binary or bounded continuous)
#'
#' @param shift_type A \code{character} defining the type of shift to be applied
#'  to the treatment. By default, this is an exponential tilt intervention.
#' @param delta A \code{numeric}, specification of the magnitude of the
#'  desired shift.
#' @param e_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a cleverly parameterized
#'  propensity score that includes the mediators, i.e., e = P(A | Z, W).
#' @param phi_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a reduced regression
#'  useful for computing the efficient one-step estimator, i.e., phi(W) =
#'  E[m(A = 1, Z, W) - m(A = 0, Z, W) | W).
#' @param max_iter A \code{numeric} setting the total number of iterations to be
#'  used in the targeted procedure based on universal least favorable submodels.
#' @param step_size A \code{numeric} giving the step size (\code{delta_epsilon}
#'  in \code{tmle3}) to be used in the targeted procedure based on universal
#'  least favorable submodels.
#' @param ... Additional arguments (currently unused).
#'
#' @export
#
tmle_medoutcon <- function(shift_type = "exptilt",
                           delta, e_learners, phi_learners,
                           max_iter = 1e4, step_size = 1e-6,
                           ...) {
  # this is a factory function
  tmle3_Spec_medoutcon$new(
    shift_type, delta,
    e_learners, phi_learners,
    max_iter, step_size,
    ...
  )
}

################################################################################

#' NPSEM for Stochastic Mediation Effects with Mediator-Outcome Confounding
#'
#' @param node_list A \code{list} object specifying the different nodes in the
#'  non-parametric structural equation model.
#' @param variable_types Used to define how variables are handled. Optional.
#'
#' @importFrom tmle3 define_node
#'
#' @keywords internal
#
medoutcon_npsem <- function(node_list, variable_types = NULL) {
  # make non-parametric structural equation model
  npsem <- list(
    # baseline covariates
    tmle3::define_node("W", node_list$W, variable_type = variable_types$W),
    # treatment
    tmle3::define_node("A", node_list$A, c("W"),
      variable_type = variable_types$A
    ),
    # mediator-outcome confounder
    tmle3::define_node("Z", node_list$Z, c("A", "W"),
      variable_type = variable_types$Z
    ),
    # mediator(s)
    tmle3::define_node("M", node_list$M, c("Z", "A", "W"),
      variable_type = variable_types$M
    ),
    # outcome
    tmle3::define_node("Y", node_list$Y, c("M", "Z", "A", "W"),
      variable_type = variable_types$Y, scale = TRUE
    )
  )
  return(npsem)
}

################################################################################

#' Stochastic Mediation Likelihood Factors
#'
#' @param tmle_task A \code{tmle3_Task} object specifying the data and the
#'  NPSEM for use in constructing elements of TML estimator.
#' @param learner_list A \code{list} specifying which learners are to be applied
#'  for each of the regression tasks required for the TML estimator.
#'
#' @importFrom tmle3 define_lf LF_emp LF_fit Likelihood
#'
#' @keywords internal
#
medoutcon_likelihood <- function(tmle_task, learner_list) {
  # covariates
  W_factor <- tmle3::define_lf(tmle3::LF_emp, "W")

  # treatment (bound likelihood away from 0 (and 1 if binary))
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continuous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }

  # treatment
  A_factor <- tmle3::define_lf(tmle3::LF_fit, "A",
    learner = learner_list[["A"]],
    bound = A_bound
  )

  # mediator-outcome confounder (bound likelihood away from 0 (and 1 if binary))
  Z_type <- tmle_task$npsem[["Z"]]$variable_type
  if (Z_type$type == "continuous") {
    Z_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (Z_type$type %in% c("binomial", "categorical")) {
    Z_bound <- 0.025
  } else {
    Z_bound <- NULL
  }

  # mediator-outcome confounder
  Z_factor <- tmle3::define_lf(tmle3::LF_fit, "Z",
    learner = learner_list[["Z"]],
    bound = Z_bound
  )

  # outcome
  Y_factor <- tmle3::define_lf(tmle3::LF_fit, "Y",
    learner = learner_list[["Y"]],
    type = "mean"
  )

  # construct and train likelihood
  factor_list <- list(W_factor, A_factor, Z_factor, Y_factor)

  likelihood_def <- tmle3::Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}

################################################################################

#' Make task for derived likelihood factor e(a | m, w)
#'
#' @param tmle_task A \code{tmle3_Task} object specifying the data and the
#'   NPSEM for use in constructing elements of TML estimator.
#' @param likelihood A trained \code{Likelihood} object from \code{tmle3},
#'  constructed via the helper function \code{stochastic_mediation_likelihood}.
#'
#' @importFrom sl3 sl3_Task
#'
#' @keywords internal
#
make_e_task <- function(tmle_task, likelihood) {
  e_data <- tmle_task$internal_data
  e_task <- sl3::sl3_Task$new(
    data = e_data,
    outcome = tmle_task$npsem[["A"]]$variables,
    covariates = c(
      tmle_task$npsem[["M"]]$variables,
      tmle_task$npsem[["W"]]$variables
    )
  )
  return(e_task)
}

################################################################################

#' Make task for derived likelihood factor r(z | a', m, w)
#'
#' @param tmle_task A \code{tmle3_Task} object specifying the data and the
#'   NPSEM for use in constructing elements of TML estimator.
#' @param likelihood A trained \code{Likelihood} object from \code{tmle3},
#'  constructed via the helper function \code{stochastic_mediation_likelihood}.
#'
#' @importFrom sl3 sl3_Task
#'
#' @keywords internal
#
make_r_task <- function(tmle_task, likelihood) {
  m_data <- tmle_task$internal_data
  m_task <- sl3::sl3_Task$new(
    data = m_data,
    outcome = tmle_task$npsem[["Z"]]$variables,
    covariates = c(
      tmle_task$npsem[["A"]]$variables,
      tmle_task$npsem[["M"]]$variables,
      tmle_task$npsem[["W"]]$variables
    )
  )
  return(r_task)
}

################################################################################

#' Make task for derived likelihood factor phi(W)
#'
#' @param tmle_task A \code{tmle3_Task} object specifying the data and the
#'  NPSEM for use in constructing elements of TML estimator.
#' @param likelihood A trained \code{Likelihood} object from \code{tmle3},
#'  constructed via the helper function \code{stochastic_mediation_likelihood}.
#'
#' @importFrom data.table as.data.table data.table
#' @importFrom uuid UUIDgenerate
#' @importFrom sl3 sl3_Task
#'
#' @keywords internal
#
make_phi_task <- function(tmle_task, likelihood) {
  # create treatment and control tasks for intervention conditions
  treatment_task <-
    tmle_task$generate_counterfactual_task(
      uuid = uuid::UUIDgenerate(),
      new_data = data.table::data.table(A = 1)
    )
  control_task <-
    tmle_task$generate_counterfactual_task(
      uuid = uuid::UUIDgenerate(),
      new_data = data.table::data.table(A = 0)
    )

  # create counterfactual outcomes and construct pseudo-outcome
  m1 <- likelihood$get_likelihood(treatment_task, "Y")
  m0 <- likelihood$get_likelihood(control_task, "Y")
  m_diff <- m1 - m0

  # create regression task for pseudo-outcome and baseline covariates
  phi_data <- data.table::as.data.table(list(
    m_diff = m_diff,
    tmle_task$get_tmle_node("W")
  ))
  phi_task <- sl3::sl3_Task$new(
    data = phi_data,
    outcome = "m_diff",
    covariates = tmle_task$npsem[["W"]]$variables,
    outcome_type = "continuous"
  )
  return(phi_task)
}
