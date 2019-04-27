#' Fit nuisance parameter u(z,a,w)
#'
#' @param data A \code{data.table} containing the observed data, with columns
#'  in the order specified by the NPSEM (Y, M, Z, A, W), with column names set
#'  appropriately based on the original input data. Such a structure is merely
#'  a convenience utility to passing data around to the various core estimation
#'  routines and is automatically generated as part of a call to the user-facing
#'  wrapper function \code{medoutcon}.
#' @param valid_data A holdout data set, with columns exactly matching those
#'  appearing in the preceding argument \code{data}, to be used for estimation
#'  via cross-fitting. Optional, defaulting to \code{NULL}.
#' @param contrast A \code{numeric} double indicating the two values of the
#'  intervention \code{A} to be compared. The default value of \code{c(0, 1)}
#'  assumes a binary intervention node \code{A}, though support for categorical
#'  interventions is planned for future releases.
#' @param lrnr_stack A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a model for the
#'  propensity score, i.e., g = P(A | W).
#' @param m_output ...
#' @param q_output ...
#' @param r_output ...
#' @param e_output ...
#'
#' @importFrom data.table as.data.table copy setnames ":="
#' @importFrom sl3 sl3_Task
#'
fit_u_param <- function(data,
                        valid_data,
                        contrast,
                        lrnr_stack,
                        m_output,
                        q_output,
                        r_output,
                        e_output) {
  # extract required elements from fit likelihood components/mechanisms
  e_natural_pred <- e_output$e_est$e_pred_A_natural
  e_star_pred <- e_output$e_est$e_pred_A_star
  q_natural_pred <- q_output$est$pred_natural
  r_natural_pred <- r_output$est$pred_natural
  m_natural_pred <- m_output$m_est$m_pred_A_natural

  # construct outcome component of nuisance parameter
  u_out <- m_natural_pred * (q_natural_pred / r_natural_pred) *
      (e_star_pred / e_natural_pred)

  # use full data for counterfactual prediction if no validation data given
  if (is.null(valid_data)) {
    # copy full data
    data_intervene <- data.table::copy(data)
  } else {
    # copy only validation data
    data_intervene <- data.table::copy(valid_data)
  }


  # construct task for nuisance parameter regression
  u_data <- data[, U := u_out]
  u_task <- sl3::sl3_Task(
    data = u_data,
    covariates = c("A", "Z", w_names),
    outcome = "U"
  )

  # fit regression model and evaluate for contrast of interest
  u_natural_fit <- lrnr_stack$train(u_task)
  u_intervened_pred <- u_natural_fit$predict(u_intervened_task)


}

################################################################################

#' Fit nuisance parameter v(a,w)
#'
#' @param data A \code{data.table} containing the observed data, with columns
#'  in the order specified by the NPSEM (Y, M, Z, A, W), with column names set
#'  appropriately based on the original input data. Such a structure is merely
#'  a convenience utility to passing data around to the various core estimation
#'  routines and is automatically generated as part of a call to the user-facing
#'  wrapper function \code{medoutcon}.
#' @param valid_data A holdout data set, with columns exactly matching those
#'  appearing in the preceding argument \code{data}, to be used for estimation
#'  via cross-fitting. Optional, defaulting to \code{NULL}.
#' @param contrast A \code{numeric} double indicating the two values of the
#'  intervention \code{A} to be compared. The default value of \code{c(0, 1)}
#'  assumes a binary intervention node \code{A}, though support for categorical
#'  interventions is planned for future releases.
#' @param lrnr_stack A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a model for the
#'  propensity score, i.e., g = P(A | W).
#' @param m_output ...
#' @param q_output ...
#'
#' @importFrom data.table as.data.table copy setnames ":="
#' @importFrom sl3 sl3_Task
#'
fit_v_param <- function(data,
                        valid_data,
                        contrast,
                        lrnr_stack,
                        m_output,
                        q_output) {
    # extract required elements from fit likelihood components/mechanisms
    q_natural_pred <- q_output$est$pred_natural
    m_natural_pred <- m_output$m_est$m_pred_A_natural

    # construct outcome component of nuisance parameter
    # NOTE: should include sum over Z
    v_out <- m_natural_pred * q_natural_pred

    # construct task for nuisance parameter regression
    v_data <- data[, V := v_out]
    v_task <- sl3::sl3_Task(
      data = v_data,
      covariates = c("A", w_names),
      outcome = "V"
    )

    # run regression...
}
