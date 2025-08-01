#' Efficient estimation of natural and interventional (in)direct effects
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar object corresponding
#'   to a set of baseline covariates.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'   parameter of interest is defined as a location shift of this quantity.
#' @param Z A \code{numeric} vector corresponding to an intermediate confounder
#'   affected by treatment (on the causal pathway between the intervention A,
#'   mediators M, and outcome Y, but unaffected itself by the mediators). When
#'   set to \code{NULL}, the natural (in)direct effects are estimated.
#' @param R A \code{logical} vector indicating whether a sampled observation's
#'   mediator was measured via a two-phase sampling design. Defaults to a
#'   vector of ones, indicating that two-phase sampling was not performed.
#' @param M A \code{numeric} vector, \code{matrix}, \code{data.frame}, or
#'   similar corresponding to a set of mediators (on the causal pathway between
#'   the intervention A and the outcome Y).
#' @param Y A \code{numeric} vector corresponding to an outcome variable.
#' @param obs_weights A \code{numeric} vector of observation-level weights. The
#'   default is to give all observations equal weighting.
#' @param svy_weights A \code{numeric} vector of observation-level weights that
#'   have been computed externally, such as survey sampling weights. Such
#'   weights are used in the construction of re-weighted efficient estimators.
#' @param two_phase_weights A \code{numeric} vector of known observation-level
#'   weights corresponding to the inverse probability of the mediator being
#'   measured. Defaults to a vector of ones.
#' @param effect A \code{character} indicating whether to compute the direct or
#'   the indirect effects of <https://doi.org/10.1093/biomet/asaa085>. This is
#'   ignored when the argument \code{contrast} is provided. By default, the
#'   direct effect is estimated.
#' @param contrast A \code{numeric} double indicating the two values of the
#'   intervention \code{A} to be compared. The default value of \code{NULL} has
#'   no effect, as the value of the argument \code{effect} is instead used to
#'   define the contrasts. To override \code{effect}, provide a \code{numeric}
#'   double vector, giving the values of a' and a*, e.g., \code{c(0, 1)}.
#' @param g_learners A \code{\link[sl3]{Stack}} object, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
#'   learners from \pkg{sl3}; used to fit a model for the propensity score.
#' @param h_learners A \code{\link[sl3]{Stack}} object, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
#'   learners from \pkg{sl3}; used to fit a model for a parameterization of the
#'   propensity score that conditions on the mediators.
#' @param b_learners A \code{\link[sl3]{Stack}} object, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
#'   learners from \pkg{sl3}; used to fit a model for the outcome regression.
#' @param q_learners A \code{\link[sl3]{Stack}} object, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
#'   learners from \pkg{sl3}; used to fit a model for a nuisance regression of
#'   the intermediate confounder, conditioning on the treatment and potential
#'   baseline covariates.
#' @param r_learners A \code{\link[sl3]{Stack}} object, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
#'   learners from \pkg{sl3}; used to fit a model for a nuisance regression of
#'   the intermediate confounder, conditioning on the mediators, the treatment,
#'   and potential baseline confounders.
#' @param u_learners A \code{\link[sl3]{Stack}} object, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
#'   learners from \pkg{sl3}; used to fit a pseudo-outcome regression required
#'   for in the efficient influence function.
#' @param v_learners A \code{\link[sl3]{Stack}} object, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
#'   learners from \pkg{sl3}; used to fit a pseudo-outcome regression required
#'   for in the efficient influence function.
#' @param d_learners A \code{\link[sl3]{Stack}} object, or other learner class
#'   (inheriting from \code{\link[sl3]{Lrnr_base}}), containing instantiated
#'   learners from \pkg{sl3}; used to fit an initial efficient influence
#'   function regression when computing the efficient influence function in a
#'   two-phase sampling design.
#' @param estimator The desired estimator of the direct or indirect effect (or
#'   contrast-specific parameter) to be computed. Both an efficient one-step
#'   estimator using cross-fitting and a cross-validated targeted minimum loss
#'   estimator (TMLE) are available. The default is the TML estimator.
#' @param estimator_args A \code{list} of additional arguments passed (via
#'   \code{...}) to the function call for the specified estimator. The default
#'   is chosen so as to allow the number of folds used to compute the one-step
#'   or TML estimators to be adjusted and for stratified cross-validation to be
#'   used in cases of rare outcomes. In the case of the TML estimator, the
#'   number of update (fluctuation) iterations is limited, and a tolerance is
#'   included for the updates introduced by tilting (fluctuation) models.
#' @param g_bounds A \code{numeric} vector containing two values, the first
#'   being the minimum allowable estimated propensity score value and the
#'   second being the maximum allowable for estimated propensity scores.
#'
#' @importFrom data.table as.data.table setnames set
#' @importFrom sl3 Lrnr_glm_fast Lrnr_hal9001
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#' # here, we show one-step and TML estimates of the interventional direct
#' # effect; the indirect effect can be evaluated by a straightforward change
#' # to the penultimate argument. the natural direct and indirect effects can
#' # be evaluated by omitting the argument Z (inappropriate in this example).
#' # create data: covariates W, exposure A, post-exposure-confounder Z,
#' #              mediator M, outcome Y
#' n_obs <- 200
#' w_1 <- rbinom(n_obs, 1, prob = 0.6)
#' w_2 <- rbinom(n_obs, 1, prob = 0.3)
#' w <- as.data.frame(cbind(w_1, w_2))
#' a <- as.numeric(rbinom(n_obs, 1, plogis(rowSums(w) - 2)))
#' z <- rbinom(n_obs, 1, plogis(rowMeans(-log(2) + w - a) + 0.2))
#' m_1 <- rbinom(n_obs, 1, plogis(rowSums(log(3) * w + a - z)))
#' m_2 <- rbinom(n_obs, 1, plogis(rowSums(w - a - z)))
#' m <- as.data.frame(cbind(m_1, m_2))
#' y <- rbinom(n_obs, 1, plogis(1 / (rowSums(w) - z + a + rowSums(m))))
#'
#' # one-step estimate of the interventional direct effect
#' os_de <- medoutcon(
#'   W = w, A = a, Z = z, M = m, Y = y,
#'   effect = "direct",
#'   estimator = "onestep"
#' )
#'
#' # TML estimate of the interventional direct effect
#' # NOTE: improved variance estimate and de-biasing from targeting procedure
#' tmle_de <- medoutcon(
#'   W = w, A = a, Z = z, M = m, Y = y,
#'   effect = "direct",
#'   estimator = "tmle"
#' )
medoutcon <- function(W,
                      A,
                      Z,
                      M,
                      Y,
                      R = rep(1, length(Y)),
                      obs_weights = rep(1, length(Y)),
                      svy_weights = NULL,
                      two_phase_weights = rep(1, length(Y)),
                      effect = c("direct", "indirect", "pm"),
                      contrast = NULL,
                      g_learners = sl3::Lrnr_glm_fast$new(),
                      h_learners = sl3::Lrnr_glm_fast$new(),
                      b_learners = sl3::Lrnr_glm_fast$new(),
                      q_learners = sl3::Lrnr_glm_fast$new(),
                      r_learners = sl3::Lrnr_glm_fast$new(),
                      u_learners = sl3::Lrnr_hal9001$new(),
                      v_learners = sl3::Lrnr_hal9001$new(),
                      d_learners = sl3::Lrnr_glm_fast$new(),
                      estimator = c("tmle", "onestep"),
                      estimator_args = list(
                        cv_folds = 10L, cv_strat = FALSE, strat_pmin = 0.1,
                        max_iter = 10L, tiltmod_tol = 5
                      ),
                      g_bounds = c(0.005, 0.995)) {
  # set defaults
  estimator <- match.arg(estimator)

  # set constant Z for estimation of the natural (in)direct effects
  if (is.null(Z)) {
    Z <- rep(1, length(Y))
    effect_type <- "natural"
  } else {
    effect_type <- "interventional"
  }
  # construct input data structure
  data <- data.table::as.data.table(cbind(
    Y, M, R, Z, A, W, obs_weights,
    two_phase_weights
  ))
  w_names <- paste("W", seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  m_names <- paste("M", seq_len(dim(data.table::as.data.table(M))[2]),
    sep = "_"
  )
  data.table::setnames(data, c(
    "Y", m_names, "R", "Z", "A", w_names,
    "obs_weights", "two_phase_weights"
  ))

  # bound outcome Y in unit interval
  min_y <- min(data[["Y"]])
  max_y <- max(data[["Y"]])
  data.table::set(data, j = "Y", value = scale_to_unit(data[["Y"]]))

  # need to loop over different contrasts to construct direct/indirect effects
  if (is.null(contrast)) {
    if (effect != "pm") {
      # select appropriate component for direct vs indirect effects
      is_effect_direct <- (effect == "direct")
      contrast_grid <- list(switch(2 - is_effect_direct,
        c(0, 0),
        c(1, 1)
      ))
      # term needed in the decomposition for both effects
      contrast_grid[[2]] <- c(1, 0)
    } else {
      contrast_grid <- list(
        c(1, 1), c(0, 0), c(1, 0)
      )
    }
  } else {
    # otherwise, simply estimate for the user-given contrast
    contrast_grid <- list(contrast)
    effect <- NULL
  }

  est_params <- lapply(contrast_grid, function(contrast) {
    if (estimator == "onestep") {
      # set arguments to pass to one-step workhorse function
      estimator_args <- estimator_args[names(estimator_args) %in%
        names(formals(est_onestep))]

      # EFFICIENT ONE-STEP ESTIMATOR
      est_onestep_args <- list(
        data = data,
        contrast = contrast,
        g_learners = g_learners,
        h_learners = h_learners,
        b_learners = b_learners,
        q_learners = q_learners,
        r_learners = r_learners,
        u_learners = u_learners,
        v_learners = v_learners,
        d_learners = d_learners,
        w_names = w_names,
        m_names = m_names,
        y_bounds = c(min_y, max_y),
        effect_type = effect_type,
        svy_weights = svy_weights,
        g_bounds = g_bounds
      )
      est_onestep_args <- unlist(list(est_onestep_args, estimator_args),
        recursive = FALSE
      )
      est_out <- do.call(est_onestep, est_onestep_args)
    } else if (estimator == "tmle") {
      # set arguments to pass to TMLE workhorse function
      estimator_args <- estimator_args[names(estimator_args) %in%
        names(formals(est_tml))]

      # TARGETED MINIMUM LOSS ESTIMATOR
      est_tml_args <- list(
        data = data,
        contrast = contrast,
        g_learners = g_learners,
        h_learners = h_learners,
        b_learners = b_learners,
        q_learners = q_learners,
        r_learners = r_learners,
        u_learners = u_learners,
        v_learners = v_learners,
        d_learners = d_learners,
        w_names = w_names,
        m_names = m_names,
        y_bounds = c(min_y, max_y),
        effect_type = effect_type,
        svy_weights = svy_weights,
        g_bounds = g_bounds
      )
      est_tml_args <- unlist(list(est_tml_args, estimator_args),
        recursive = FALSE
      )
      est_out <- do.call(est_tml, est_tml_args)
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
      type = estimator,
      param = paste("direct", effect_type, sep = "_"),
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
      type = estimator,
      param = paste("indirect", effect_type, sep = "_"),
      outcome = as.numeric(Y)
    )
    class(ie_est_out) <- "medoutcon"
    return(ie_est_out)
  } else if (is.null(contrast) && (effect == "pm")) {
    # compute parameter estimate, influence function, and variances
    pm_theta_est <- 1 - log(est_params[[3]]$theta / est_params[[2]]$theta) /
      log(est_params[[1]]$theta / est_params[[2]]$theta)
    pm_eif_est <- -est_params[[3]]$eif /
      (est_params[[3]]$theta_plugin * log(est_params[[1]]$theta_plugin /
        est_params[[2]]$theta_plugin)) +
      est_params[[2]]$eif * (
        (log(est_params[[1]]$theta_plugin / est_params[[2]]$theta_plugin) -
          log(est_params[[3]]$theta_plugin / est_params[[2]]$theta_plugin)) /
          (est_params[[2]]$theta_plugin *
            (log(est_params[[1]]$theta_plugin /
              est_params[[2]]$theta_plugin))^2)) +
      est_params[[1]]$eif * log(est_params[[3]]$theta_plugin /
        est_params[[2]]$theta_plugin) /
        (est_params[[1]]$theta_plugin * (log(est_params[[1]]$theta_plugin /
          est_params[[2]]$theta_plugin))^2)
    pm_var_est <- stats::var(pm_eif_est) / nrow(data)

    # construct output in same style as for contrast-specific parameter
    pm_est_out <- list(
      theta = pm_theta_est,
      var = pm_var_est,
      eif = pm_eif_est,
      type = estimator,
      param = paste("pm", effect_type, sep = "_"),
      outcome = as.numeric(Y),
      contrast_results = list(
        contrast_1_1_mean = est_params[[1]]$theta,
        contrast_1_1_eif = est_params[[1]]$eif,
        contrast_0_0_mean = est_params[[2]]$theta,
        contrast_0_0_eif = est_params[[2]]$eif,
        contrast_1_0_mean = est_params[[3]]$theta,
        contrast_1_0_eif = est_params[[3]]$eif
      )
    )
    class(pm_est_out) <- "medoutcon"
    return(pm_est_out)
  } else {
    est_out <- unlist(est_params, recursive = FALSE)
    est_out$param <- paste0(
      "tsm(", "a'=", contrast[1], ",",
      "a*=", contrast[2], ")"
    )
    est_out$.contrast <- contrast
    class(est_out) <- "medoutcon"
    return(est_out)
  }
}
