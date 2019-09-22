#' Parameter for mediation effect under stochastic interventions
#'
#' Parameter definition class. See https://arxiv.org/abs/1901.02776
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom tmle3 Param_base
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_medoutcon, shift_param, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}}
#'           corresponding to the observed likelihood.
#'     }
#'     \item{\code{shift_param}}{\code{numeric}, specification of the magnitude
#'           of the desired shift (a multiplier for the propensity score).
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be
#'           treated as the outcome
#'     }
#'   }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood}}{the counterfactual likelihood under the
#'           stochastic intervention on mediators and treatment.
#'     }
#'     \item{\code{lf_exptilt}}{Object derived from \code{\link{LF_base}} for
#'           assessing the stochastic mediation intervention.
#'     }
#'     \item{\code{treatment_task}}{\code{\link{tmle3_Task}} object created from
#'           setting the intervention to the treatment condition: do(A = 1).
#'     }
#'     \item{\code{control_task}}{\code{\link{tmle3_Task}} object created from
#'           setting the intervention to the control condition: do(A = 0).
#'     }
#'     \item{\code{shift_param}}{\code{numeric}, specification of the magnitude
#'           of the desired shift (a multiplier for the propensity score).
#'     }
#' }
#' @export
Param_medshift <- R6::R6Class(
  classname = "Param_medoutcon",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3::Param_base,
  public = list(
    initialize = function(observed_likelihood,
                              shift_param,
                              ...,
                              outcome_node = "Y") {
      # copied from standard parameter definitions
      super$initialize(observed_likelihood, list(...),
        outcome_node = outcome_node
      )
      tmle_task <- observed_likelihood$training_task
      # counterfactual tasks
      treatment_task <-
        tmle_task$generate_counterfactual_task(
          uuid = uuid::UUIDgenerate(),
          new_data = data.table(A = 1)
        )
      control_task <-
        tmle_task$generate_counterfactual_task(
          uuid = uuid::UUIDgenerate(),
          new_data = data.table(A = 0)
        )

      # generate counterfactual likelihood under intervention via LF_exptilt
      lf_exptilt <- LF_exptilt_ipsi$new(
        name = "A",
        likelihood_base = observed_likelihood,
        shift_param = shift_param,
        treatment_task = treatment_task,
        control_task = control_task,
        cache = FALSE
      )

      # store components
      private$.cf_likelihood <- CF_Likelihood$new(
        observed_likelihood,
        lf_exptilt
      )
      private$.lf_exptilt <- lf_exptilt
      private$.shift_param <- shift_param
      private$.treatment_task <- treatment_task
      private$.control_task <- control_task
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      # get observed likelihood
      likelihood <- self$observed_likelihood
      cf_likelihood <- self$cf_likelihood
      shift_param <- self$shift_param
      treatment_task <- self$treatment_task
      control_task <- self$control_task

      # extract various likelihood components
      m_est <- likelihood$get_likelihood(tmle_task, "Y", fold_number)
      e_est <- likelihood$get_likelihood(tmle_task, "E", fold_number)
      phi_est <- likelihood$get_likelihood(tmle_task, "phi", fold_number)
      g_delta_est <- cf_likelihood$get_likelihood(tmle_task, "A", fold_number)

      # compute/extract g(1|W) for clever covariate for score of A
      g1_est <- likelihood$get_likelihood(treatment_task, "A", fold_number)
      g0_est <- likelihood$get_likelihood(control_task, "A", fold_number)

      # clever covariates
      HY <- g_delta_est / e_est

      # NOTE: exp(shift_param) for generalized exponential tilting
      HA <- (shift_param * phi_est) / ((shift_param * g1_est) + g0_est)^2

      # output clever covariates
      return(list(Y = HY, A = HA))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      # get observed likelihood
      likelihood <- self$observed_likelihood
      cf_likelihood <- self$cf_likelihood
      shift_param <- self$shift_param
      treatment_task <- self$treatment_task
      control_task <- self$control_task

      # extract various likelihood components
      y <- tmle_task$get_tmle_node(self$outcome_node)
      a <- tmle_task$get_tmle_node(self$lf_exptilt$name)
      m_est <- likelihood$get_likelihood(tmle_task, "Y", fold_number)

      # compute/extract g(1|W) for clever covariate for score of A
      g1_est <- likelihood$get_likelihood(treatment_task, "A", fold_number)
      g1_delta_est <- cf_likelihood$get_likelihood(
        treatment_task, "A",
        fold_number
      )
      g0_delta_est <- cf_likelihood$get_likelihood(
        control_task, "A",
        fold_number
      )
      m1_est <- likelihood$get_likelihood(treatment_task, "Y", fold_number)
      m0_est <- likelihood$get_likelihood(control_task, "Y", fold_number)

      # clever_covariates happens here but this is repeated computation
      HY <- self$clever_covariates(
        tmle_task,
        fold_number
      )[[self$outcome_node]]
      HA <- self$clever_covariates(
        tmle_task,
        fold_number
      )[[self$lf_exptilt$name]]

      # compute individual scores for DY, DA, DZW
      D_Y <- HY * (y - m_est)
      D_A <- HA * (a - g1_est)
      D_ZW <- (g1_delta_est * m1_est) + (g0_delta_est * m0_est)

      # parameter and influence function
      theta <- mean(D_ZW)
      eif <- D_Y + D_A + D_ZW - theta

      # output
      result <- list(psi = theta, IC = eif)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf(
        "E[%s_{%s}]", self$outcome_node,
        self$cf_likelihood$name
      )
      return(param_form)
    },
    cf_likelihood = function() {
      return(private$.cf_likelihood)
    },
    lf_exptilt = function() {
      return(private$.lf_exptilt)
    },
    shift_param = function() {
      return(private$.shift_param)
    },
    treatment_task = function() {
      return(private$.treatment_task)
    },
    control_task = function() {
      return(private$.control_task)
    },
    update_nodes = function() {
      # TODO: stop hard-coding A everywhere
      return(c(self$outcome_node, "A"))
    }
  ),
  private = list(
    .type = "medoutcon",
    .cf_likelihood = NULL,
    .lf_exptilt = NULL,
    .shift_param = NULL,
    .treatment_task = NULL,
    .control_task = NULL
  )
)
