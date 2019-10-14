#' Efficient estimation of stochastic (in)direct effects in the presence of
#' mediator-outcome confounders affected by exposure
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar object corresponding
#'  to a set of baseline covariates.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param Z A \code{numeric} vector corresponding to a mediator-outcome
#'  confounder affected by treatment (on the causal pathway between intervention
#'  A, mediator M, and outcome Y, but unaffected itself by the mediator M).
#' @param M A \code{numeric} vector, \code{matrix}, \code{data.frame}, or
#'  similar corresponding to a set of mediators (on the causal pathway between
#'  the intervention A and the outcome Y).
#' @param Y A \code{numeric} vector corresponding to an outcome variable.
#' @param obs_weights A \code{numeric} vector of observation-level weights to be
#'  incorporated in all procedures estimating nuisance parameters. The default
#'  is to give all observations equal weighting.
#' @param ext_weights A \code{numeric} vector of observation-level weights that
#'  have been computed externally. Such weights are used in the construction of
#'  a re-weighted one-step estimator or in solving a re-weighted estimating
#'  equation in the case of the TML estimator. Note that, unlike the argument
#'  \code{obs_weights}, these weights are not incorporated in the estimation of
#'  nuisance parameters. Input weights should be normalized. Use with caution.
#' @param effect A \code{character} indicating whether to compute the direct
#'  effect or the indirect effect as discussed in CITE PAPER. Note that this is
#'  ignored when the argument \code{contrast} is provided.
#' @param contrast A \code{numeric} double indicating the two values of the
#'  intervention \code{A} to be compared. The default value of \code{NULL} has
#'  no effect, as the value of the argument \code{effect} is instead used to
#'  define the contrasts. To override \code{effect}, provide a \code{numeric}
#'  double vector, giving the values of a' and a* (e.g., \code{c(0, 1)}.
#' @param g_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, for use in fitting a model for the propensity
#'  score, i.e., \eqn{g = P(A | W)}.
#' @param e_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, for use in fitting a cleverly parameterized
#'  propensity score that includes the mediators, i.e., \eqn{e = P(A | M, W)}.
#' @param m_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting the outcome regression.
#' @param q_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a regression involving
#'  the mediator-outcome confounder, i.e., \eqn{q(Z | A, W)}.
#' @param r_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a regression involving
#'  the mediator-outcome confounder, i.e., \eqn{r(Z | A, M, W)}.
#' @param u_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a reduced regression
#'  useful for computing the efficient one-step estimator, i.e.,
#'  \eqn{u(Z, A, W) = E[m(A, Z, M, W) * (q(Z | A, W) / r(Z | A, M, W)) *
#'  (e(a' | M, W) / e(A | M, W)) | Z = z, A = a, W = w]}.
#' @param v_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a reduced regression
#'  useful for computing the efficient one-step estimator, i.e.,
#'  \eqn{v(A, W) = E[\int_z m(a', z, M, W) * q(z | a', W) d \nu(z) |
#'  A = a, W = w)]}.
#' @param estimator The desired estimator of the direct or indirect effect (or
#'   contrast-specific parameter) to be computed. Currently, the singular option
#'   is an efficient one-step estimator, though support is planned for classical
#'   estimators (substitution; inverse probability weighted) and an efficient
#'   targeted maximum likelihood estimator.
#' @param estimator_args A \code{list} of extra arguments to be passed (via
#'  \code{...}) to the function call for the specified estimator. The default
#'  is so chosen as to allow the number of folds used in computing the one-step
#'  estimator to be easily adjusted.
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom stats binomial var
#' @importFrom sl3 Lrnr_glm_fast
#'
#' @export
#
medoutcon <- function(W,
                      A,
                      Z,
                      M,
                      Y,
                      obs_weights = rep(1, length(Y)),
                      ext_weights = NULL,
                      effect = c("direct", "indirect"),
                      contrast = NULL,
                      g_learners =
                        sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                      e_learners =
                        sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                      m_learners = sl3::Lrnr_glm_fast$new(),
                      q_learners =
                        sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                      r_learners =
                        sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                      u_learners = sl3::Lrnr_glm_fast$new(),
                      v_learners = sl3::Lrnr_glm_fast$new(),
                      estimator = c(
                        "onestep",
                        "tmle",
                        "ipw",
                        "sub"
                      ),
                      estimator_args = list(cv_folds = 10)) {
  # set defaults
  effect <- match.arg(effect)
  estimator <- match.arg(estimator)
  estimator_args <- unlist(estimator_args, recursive = FALSE)

  # construct input data structure
  data <- data.table::as.data.table(cbind(Y, M, Z, A, W, obs_weights))
  w_names <- paste("W", seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  m_names <- paste("M", seq_len(dim(data.table::as.data.table(M))[2]),
    sep = "_"
  )
  data.table::setnames(data, c("Y", m_names, "Z", "A", w_names, "obs_weights"))

  # need to loop over different contrasts to construct direct/indirect effects
  if (is.null(contrast)) {
    # select appropriate component for direct vs indirect effects
    is_effect_direct <- (effect == "direct")
    contrast_grid <- list(switch(2 - is_effect_direct, c(0, 0), c(1, 1)))
    # this term is needed in the decomposition for both effects
    contrast_grid[[2]] <- c(1, 0)
  } else {
    # otherwise, simply estimate for the user-given contrast
    contrast_grid <- list(contrast)
    effect <- NULL
  }

  est_params <- lapply(contrast_grid, function(contrast) {
    if (estimator == "sub") {
      # SUBSTITUTION ESTIMATOR
      stop("The substitution estimator is currently under development.")
    } else if (estimator == "ipw") {
      # INVERSE PROBABILITY RE-WEIGHTED ESTIMATOR
      stop("The IPW estimator is currently under development.")
    } else if (estimator == "onestep") {
      # EFFICIENT ONE-STEP ESTIMATOR
      onestep_est_args <- list(
        data = data,
        contrast = contrast,
        g_learners = g_learners,
        e_learners = e_learners,
        m_learners = m_learners,
        q_learners = q_learners,
        r_learners = r_learners,
        u_learners = u_learners,
        v_learners = v_learners,
        w_names = w_names,
        m_names = m_names,
        ext_weights = ext_weights,
        estimator_args
      )
      est_out <- do.call(est_onestep, onestep_est_args)
    } else if (estimator == "tmle") {
      # TARGETED MINIMUM LOSS ESTIMATOR
      tmle_est_args <- list(
        data = data,
        contrast = contrast,
        g_learners = g_learners,
        e_learners = e_learners,
        m_learners = m_learners,
        q_learners = q_learners,
        r_learners = r_learners,
        u_learners = u_learners,
        v_learners = v_learners,
        w_names = w_names,
        m_names = m_names,
        estimator_args
      )
      est_out <- do.call(est_tml, tmle_est_args)
    }

    # lazily create output as classed list
    est_out$outcome <- as.numeric(Y)
    class(est_out) <- "medoutcon"
    return(est_out)
  })

  # put effects together
  if (is.null(contrast) && (effect == "direct")) {
    # compute parameter estimate, influence function, and variances
    de_theta_est <- est_params[[2]]$theta - est_params[[1]]$theta
    de_eif_est <- est_params[[2]]$eif - est_params[[1]]$eif
    de_var_est <- stats::var(de_eif_est) / nrow(data)

    # construct output in same style as for contrast-specific parameter
    de_est_out <- list(
      theta = de_theta_est,
      var = de_var_est,
      eif = de_eif_est,
      type = "onestep",
      param = "direct_effect",
      outcome = as.numeric(Y)
    )
    class(de_est_out) <- "medoutcon"
    return(de_est_out)
  } else if (is.null(contrast) && (effect == "indirect")) {
    # compute parameter estimate, influence function, and variances
    ie_theta_est <- est_params[[1]]$theta - est_params[[2]]$theta
    ie_eif_est <- est_params[[1]]$eif - est_params[[2]]$eif
    ie_var_est <- stats::var(ie_eif_est) / nrow(data)

    # construct output in same style as for contrast-specific parameter
    ie_est_out <- list(
      theta = ie_theta_est,
      var = ie_var_est,
      eif = ie_eif_est,
      type = "onestep",
      param = "indirect_effect",
      outcome = as.numeric(Y)
    )
    class(ie_est_out) <- "medoutcon"
    return(ie_est_out)
  } else {
    est_out <- unlist(est_params, recursive = FALSE)
    est_out$param <- "contrast_spec"
    class(est_out) <- "medoutcon"
    return(est_out)
  }
}
