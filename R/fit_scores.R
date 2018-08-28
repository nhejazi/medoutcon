utils::globalVariables(c("..w_names"))
################################################################################

#' Get Dzw component of efficient influence function from nuisance parameters
#'
#' @param g_output Object containing results from fitting the propensity score
#'  regression, as produced by a call to \code{fit_g_mech}.
#' @param m_output Object containing results from fitting the outcome
#'  regression, as produced by a call to \code{fit_m_mech}.
#' @param shift_type A choice of the type of stochastic treatment regime to use
#'  -- either \code{"mtp"} for a modified treatment policy that shifts the
#'  center of the observed intervention distribution by the scalar \code{delta}
#'  or \code{"ipsi"} for an incremental propensity score shift that multiples
#'  the odds of receiving the intervention by the scalar \code{delta}.
#'
#' @keywords internal
#
compute_Dzw <- function(g_output,
                        m_output,
                        shift_type = c("ipsi", "mtp")) {
  # set IPSI shift as default for now...
  shift_type <- match.arg(shift_type)

  if (shift_type == "ipsi") {
    # get g components from output for that nuisance parameter
    g_shifted_A1 <- g_output$g_est$g_pred_shifted_A1
    g_shifted_A0 <- g_output$g_est$g_pred_shifted_A0

    # get m components from output for that nuisance parameter
    m_pred_A1 <- m_output$m_pred$m_pred_A1
    m_pred_A0 <- m_output$m_pred$m_pred_A0

    # compute component Dzw from nuisance parameters
    Dzw_A1 <- g_shifted_A1 * m_pred_A1
    Dzw_A0 <- g_shifted_A0 * m_pred_A0

    # output as simple list
    return(list(
      dzw_cntrl = Dzw_A0,
      dzw_treat = Dzw_A1
    ))
  } else if (shift_type == "mtp") {
    # approximate Monte Carlo integral using inverse uniform weighting
    Dzw_int <- integrate_over_g(g_mech = g_output$g_est$g_shifted,
                                a_vals = g_output$a_vals$a_shifted,
                                weighting = m_output$m_pred$m_natural)

    # output as simple list
    return(list(
      dzw = Dzw_int
    ))
  }
}

################################################################################

#' Get inverse probability weighted (IPW) estimate from nuisance parameters
#'
#' @param g_output Object containing results from fitting the propensity score
#'  regression, as produced by a call to \code{fit_g_mech}.
#' @param e_output Object containing results from fitting the propensity score
#'  regression while conditioning on mediators, as produced by a call to
#'  \code{fit_e_mech}.
#' @param idx_treat A \code{numeric} vector providing the indices corresponding
#'  to units that received treatment (A = 1), for a binary intervention.
#' @param idx_cntrl A \code{numeric} vector providing the indices corresponding
#'  to units that did not receive treatment (A = 0), for a binary intervention.
#' @param shift_type A choice of the type of stochastic treatment regime to use
#'  -- either \code{"mtp"} for a modified treatment policy that shifts the
#'  center of the observed intervention distribution by the scalar \code{delta}
#'  or \code{"ipsi"} for an incremental propensity score shift that multiples
#'  the odds of receiving the intervention by the scalar \code{delta}.
#'
#' @keywords internal
#
compute_ipw <- function(g_output,
                        e_output,
                        idx_treat = NULL,
                        idx_cntrl = NULL,
                        shift_type = c("ipsi", "mtp")) {
  # set IPSI shift as default for now...
  shift_type <- match.arg(shift_type)

  if (shift_type == "ipsi") {
    # extract components for A = 0 and A = 1 cases
    g_shifted_A1 <- g_output$g_est$g_pred_shifted_A1
    g_shifted_A0 <- g_output$g_est$g_pred_shifted_A0
    e_pred_A1 <- e_output$e_est$e_pred_A1
    e_pred_A0 <- e_output$e_est$e_pred_A0

    # subset computed components based on observed treatment status for g
    g_shifted_obs <- rep(NA, length(idx_treat) + length(idx_cntrl))
    g_shifted_A1_obs <- g_shifted_A1[idx_treat]
    g_shifted_A0_obs <- g_shifted_A0[idx_cntrl]
    g_shifted_obs[idx_treat] <- g_shifted_A1_obs
    g_shifted_obs[idx_cntrl] <- g_shifted_A0_obs

    # subset computed components based on observed treatment status for e
    e_pred_obs <- rep(NA, length(idx_treat) + length(idx_cntrl))
    e_pred_A1_obs <- e_pred_A1[idx_treat]
    e_pred_A0_obs <- e_pred_A0[idx_cntrl]
    e_pred_obs[idx_treat] <- e_pred_A1_obs
    e_pred_obs[idx_cntrl] <- e_pred_A0_obs

    # Hajek/stabilized weights by dividing by sample average since E[g/e] = 1
    mean_ipw <- mean(g_shifted_obs / e_pred_obs)

    # output as simple list
    return(list(
      g_shifted = g_shifted_obs,
      e_pred = e_pred_obs,
      mean_ipw = mean_ipw
    ))
  } else if (shift_type == "mtp") {
    # extract components of g and e mechanisms based on intervention type
    g_shifted_pred <- g_output$g_est$g_shifted
    e_natural_pred <- e_output$e_est$e_natural

    # Hajek/stabilized weights by dividing by sample average since E[g/e] = 1
    mean_ipw <- mean(g_shifted_pred / e_natural_pred)

    # output as simple list
    return(list(
      g_shifted = g_shifted_pred,
      e_pred = e_natural_pred,
      mean_ipw = mean_ipw
    ))
  }
}

################################################################################

#' Cross-validated evaluation of EIF components
#'
#' @param fold Object specifying cross-validation folds as generated by a call
#'  to \code{origami::make_folds}.
#' @param data A \code{data.table} containing the observed data, with columns
#'  in the order specified by the NPSEM (Y, Z, A, W), with column names set
#'  appropriately based on the original input data. Such a structure is merely
#'  a convenience utility to passing data around to the various core estimation
#'  routines and is automatically generated as part of a call to the user-facing
#'  wrapper function \code{medshift}.
#' @param delta A \code{numeric} value indicating the degree of shift in the
#'  intervention to be used in defining the causal quantity of interest. In the
#'  case of binary interventions, this takes the form of an incremental
#'  propensity score shift, acting as a multiplier of the probability with which
#'  a given observational unit receives the intervention (EH Kennedy, 2018,
#'  JASA; <doi:10.1080/01621459.2017.1422737>).
#' @param lrnr_stack_g A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a model for the
#'  propensity score, i.e., g = P(A | W).
#' @param lrnr_stack_e A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a cleverly parameterized
#'  propensity score that includes the mediators, i.e., e = P(A | Z, W).
#' @param lrnr_stack_m A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting the outcome regression,
#'  i.e., m(A, Z, W).
#' @param lrnr_stack_phi A \code{Stack} object, or other learner class
#'  (inheriting from \code{Lrnr_base}), containing a single or set of
#'  instantiated learners from the \code{sl3} package, to be used in fitting a
#'  reduced regression useful for computing the efficient one-step estimator,
#'  i.e., phi(W) = E[m(A = 1, Z, W) - m(A = 0, Z, W) | W).
#' @param w_names A \code{character} vector of the names of the columns that
#'  correspond to baseline covariates (W). The input for this argument is
#'  automatically generated by a call to the wrapper function \code{medshift}.
#' @param z_names A \code{character} vector of the names of the columns that
#'  correspond to mediators (Z). The input for this argument is automatically
#'  generated by a call to the wrapper function \code{medshift}.
#' @param shift_type A choice of the type of stochastic treatment regime to use
#'  -- either \code{"mtp"} for a modified treatment policy that shifts the
#'  center of the observed intervention distribution by the scalar \code{delta}
#'  or \code{"ipsi"} for an incremental propensity score shift that multiples
#'  the odds of receiving the intervention by the scalar \code{delta}.
#'
#' @importFrom data.table data.table
#' @importFrom origami training validation fold_index
#'
#' @keywords internal
#
cv_eif <- function(fold,
                   data,
                   delta,
                   lrnr_stack_g,
                   lrnr_stack_e,
                   lrnr_stack_m,
                   lrnr_stack_phi,
                   w_names,
                   z_names,
                   shift_type = c("ipsi", "mtp")) {
  # set IPSI shift as default for now...
  shift_type <- match.arg(shift_type)

  # make training and validation data
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  # 1) fit regression for incremental propensity score intervention
  g_out <- fit_g_mech(
    data = train_data, valid_data = valid_data,
    delta = delta,
    lrnr_stack = lrnr_stack_g, w_names = w_names,
    shift_type = shift_type
  )

  # 2) fit clever regression for treatment, conditional on mediators
  e_out <- fit_e_mech(
    data = train_data, valid_data = valid_data,
    lrnr_stack = lrnr_stack_e,
    z_names = z_names, w_names = w_names,
    shift_type = shift_type
  )

  # 3) fit regression for incremental propensity score intervention
  m_out <- fit_m_mech(
    data = train_data, valid_data = valid_data,
    lrnr_stack = lrnr_stack_m,
    z_names = z_names, w_names = w_names,
    shift_type = shift_type
  )

  # 4) difference-reduced dimension regression for phi
  if (shift_type == "ipsi") {
    phi_est <- fit_phi_mech_ipsi(
      data = valid_data, lrnr_stack = lrnr_stack_phi,
      m_output = m_out, w_names = w_names
    )
  } else if (shift_type == "mtp") {
    phi_est <- fit_phi_mech_mtp(
      data = valid_data, lrnr_stack = lrnr_stack_phi,
      m_output = m_out, e_output = e_out,
      g_output = g_out, w_names = w_names
    )
  }

  if (shift_type == "ipsi") {
    # get indices of treated and control units in validation data
    idx_A1 <- which(valid_data$A == 1)
    idx_A0 <- which(valid_data$A == 0)

    # compute component Dzw from nuisance parameters
    Dzw_groupwise <- compute_Dzw(g_output = g_out, m_output = m_out,
                                 shift_type = shift_type)
    Dzw <- Dzw_groupwise$dzw_cntrl + Dzw_groupwise$dzw_treat

    # compute component Da from nuisance parameters
    g_pred_A1 <- g_out$g_est$g_pred_A1
    g_pred_A0 <- g_out$g_est$g_pred_A0
    Da_numerator <- delta * phi_est * (valid_data$A - g_pred_A1)
    Da_denominator <- (delta * g_pred_A1 + g_pred_A0)^2
    Da <- Da_numerator / Da_denominator

    # compute component Dy from nuisance parameters
    ipw_groupwise <- compute_ipw(
      g_output = g_out, e_output = e_out,
      idx_treat = idx_A1, idx_cntrl = idx_A0,
      shift_type = shift_type
    )

    # stabilize weights in AIPW by dividing by sample average since E[g/e] = 1
    mean_ipw <- ipw_groupwise$mean_ipw
    g_shifted <- ipw_groupwise$g_shifted
    e_pred <- ipw_groupwise$e_pred
    sipw <- ((g_shifted / e_pred) / mean_ipw)

    # extract outcome component mechanism for estimting Dy
    m_pred_obs <- rep(NA, nrow(valid_data))
    m_pred_A1_obs <- m_out$m_pred$m_pred_A1[idx_A1]
    m_pred_A0_obs <- m_out$m_pred$m_pred_A0[idx_A0]
    m_pred_obs[idx_A1] <- m_pred_A1_obs
    m_pred_obs[idx_A0] <- m_pred_A0_obs

    # assemble Dy component estimate
    Dy <- sipw * (valid_data$Y - m_pred_obs)

  } else if (shift_type == "mtp") {
    # compute component Dzw from nuisance parameters
    Dzw_est <- compute_Dzw(g_output = g_out, m_output = m_out,
                           shift_type = shift_type)
    Dzw <- rep(as.numeric(sum(Dzw_est$dzw)), times = nrow(g_out$g_est))

    # compute component Da from nuisance parameters
    g_natural <- g_out$g_est$g_natural
    a_natural <- g_out$a_vals$a_natural

    # approximate Monte Carlo integral using inverse uniform weighting
    int_Da_phi <- integrate_over_g(g_mech = g_natural,
                                   a_vals = a_natural,
                                   weighting = phi_est)
    Da <- phi_est - sum(int_Da_phi)

    # compute component Dy from nuisance parameters
    ipw_out <- compute_ipw(
      g_output = g_out, e_output = e_out,
      shift_type = shift_type
    )

    # stabilize weights in AIPW by dividing by sample average since E[g/e] = 1
    mean_ipw <- ipw_out$mean_ipw
    g_shifted <- ipw_out$g_shifted
    e_pred <- ipw_out$e_pred
    sipw <- ((g_shifted / e_pred) / mean_ipw)

    # extract outcome mechanism estimate under natural intervention value
    m_pred <- m_out$m_pred$m_natural

    # assemble Dy component estimate
    Dy <- sipw * (valid_data$Y - m_pred)
  }

  # output list
  out <- list(data.table::data.table(
    Dy = Dy, Da = Da, Dzw = Dzw,
    fold = origami::fold_index()
  ))
  return(out)
}
