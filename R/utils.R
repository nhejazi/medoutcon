utils::globalVariables(c("effect", "param"))

#' Confidence intervals for interventional mediation effect estimates
#'
#' Compute confidence intervals for objects of class \code{medoutcon}, which
#' contain estimates produced by \code{\link{medoutcon}}.
#'
#' @param object An object of class \code{medoutcon}, as produced by invoking
#'  \code{\link{medoutcon}}, for which a confidence interval is to be computed.
#' @param parm A \code{numeric} vector indicating indices of \code{object$est}
#'  for which to return confidence intervals.
#' @param level A \code{numeric} indicating the level of the confidence
#'  interval to be computed.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats qnorm qlogis plogis
#' @importFrom assertthat assert_that
#'
#' @method confint medoutcon
#'
#' @export
confint.medoutcon <- function(object,
                              parm = seq_len(object$theta),
                              level = 0.95,
                              ...) {
  # inference is currently limited to the efficient estimators
  assertthat::assert_that(object$type %in% c("onestep", "tmle"))

  # first, let's get Z_(1 - alpha)
  ci_norm_bounds <- c(-1, 1) * abs(stats::qnorm(p = (1 - level) / 2))

  # assume continuous outcome if more than two levels in outcome node
  if (length(unique(object$outcome)) > 2 ||
    stringr::str_detect(object$param, "direct")) {
    # NOTE: variance already scaled (i.e., Var(D)/n)
    se_eif <- sqrt(object$var)

    # compute the interval around the point estimate
    ci_theta <- ci_norm_bounds * se_eif + object$theta
  } else if (length(unique(object$outcome)) == 2) {
    # for binary outcomes, create CI on the logit scale and back-transform
    theta_ratio <- stats::qlogis(object$theta)
    grad_ratio_delta <- (1 / object$theta) + (1 / (1 - object$theta))
    se_eif_logit <- sqrt(grad_ratio_delta^2 * object$var)
    ci_theta <- stats::plogis(ci_norm_bounds * se_eif_logit + theta_ratio)
  }

  # set up output CI object
  ci_out <- c(ci_theta[1], object$theta, ci_theta[2])
  names(ci_out) <- c("lwr_ci", "est", "upr_ci")
  return(ci_out)
}

###############################################################################

#' Summary for interventional mediation effect estimate objects
#'
#' Print a convenient summary for objects of \code{S3} class \code{medoutcon}.
#'
#' @param object An object of class \code{medoutcon}, as produced by invoking
#'  \code{\link{medoutcon}}.
#' @param ... Other arguments. Not currently used.
#' @param ci_level A \code{numeric} indicating the level of the confidence
#'  interval to be computed.
#'
#' @method summary medoutcon
#'
#' @importFrom stats confint
#' @importFrom tibble as_tibble
#'
#' @export
summary.medoutcon <- function(object,
                              ...,
                              ci_level = 0.95) {
  # compute confidence interval
  est_with_ci <- stats::confint(object, level = ci_level)

  # create output table from input object and confidence interval results
  est_summary <- tibble::as_tibble(list(
    lwr_ci = est_with_ci[1],
    param_est = est_with_ci[2],
    upr_ci = est_with_ci[3],
    var_est = object$var,
    eif_mean = mean(object$eif),
    estimator = object$type,
    param = object$param
  ))

  # store CI level as hidden attribute
  attr(est_summary, ".ci_level") <- ci_level
  return(est_summary)
}

###############################################################################

#' Print method for interventional mediation effect estimate objects
#'
#' The \code{print} method for objects of class \code{medoutcon}.
#'
#' @param x An object of class \code{medoutcon}.
#' @param ... Other options (not currently used).
#'
#' @method print medoutcon
#'
#' @importFrom zeallot "%<-%"
#' @importFrom scales percent
#' @importFrom stringr str_detect str_split str_to_title
#'
#' @export
print.medoutcon <- function(x, ...) {
  # get summary, including confidence interval
  x_summary <- summary(x)
  ci_level <- attr(x_summary, ".ci_level")

  # construct and print output
  if (stringr::str_detect(x$param, "tsm")) {
    # TODO: printing specific to counterfactual mean
    message("Counterfactual TSM")
    message(paste0(
      "Contrast: A = ", x$.contrast[1], ", ",
      paste0("M(A = ", x$.contrast[2], ")")
    ))
  } else {
    # mangle the abbreviated parameter name
    c(param, effect) %<-% unlist(stringr::str_split(x_summary$param, "_"))
    message(stringr::str_to_title(paste(effect, param, "effect")))
  }
  message("Estimator: ", x_summary$estimator)
  message("Estimate: ", round(x_summary$param_est, 3))
  message("Std. Error: ", round(sqrt(x_summary$var_est), 3))
  message(
    scales::percent(ci_level), " CI: [",
    round(x_summary$lwr_ci, 3), ", ", round(x_summary$upr_ci, 3), "]"
  )
}

###############################################################################

#' Bounding to numerical precision
#'
#' Bounds extreme values to a specified tolerance level, for use with sensitive
#' quantities that must be transformed, e.g., via \code{\link[stats]{qlogis}}.
#'
#' @param vals A \code{numeric} vector of values in the interval [0, 1].
#' @param tol A \code{numeric} indicating the tolerance limit to which extreme
#'  values should be truncated. Realizations of \code{val} less than \code{tol}
#'  are truncated to \code{tol} while those greater than (1 - \code{tol}) are
#'  truncated to (1 - \code{tol}).
#'
#' @importFrom assertthat assert_that
#'
#' @keywords internal
bound_precision <- function(vals, tol = 1e-6) {
  vals[vals < tol] <- tol
  vals[vals > 1 - tol] <- 1 - tol
  return(vals)
}

###############################################################################

#' Bounding propensity scores
#'
#' Bounds estimated propensity score values to be within a specified range.
#'
#' @param vals A \code{numeric} vector of values in the interval [0, 1].
#' @param bounds A \code{numeric} vector containing two values, the first being
#'  the minimum allowable value and the second being the maximum allowable for
#'  values appearing in the vector \code{vals} (the previous argument).
#'
#' @importFrom assertthat assert_that
#'
#' @keywords internal
bound_propensity <- function(vals, bounds = c(0.001, 0.999)) {
  assertthat::assert_that(!(max(vals) > 1 || min(vals) < 0))
  vals[vals < bounds[1]] <- bounds[1]
  vals[vals > bounds[2]] <- bounds[2]
  return(vals)
}

###############################################################################

#' Scale values to the unit interval
#'
#' @param vals A \code{numeric} vector of values to be scaled into the closed
#'  interval [0, 1].
#'
#' @keywords internal
scale_to_unit <- function(vals) {
  vals_scaled <- (vals - min(vals)) / (max(vals) - min(vals))
  return(vals_scaled)
}

###############################################################################

#' Scale values from the unit interval to their original scale
#'
#' @param scaled_vals A \code{numeric} vector of values scaled to lie in the
#'  closed interval [0, 1] by use of \code{\link{scale_to_unit}}.
#' @param max_orig The maximum of the values on the original scale.
#' @param min_orig The minimum of the values on the original scale.
#'
#' @keywords internal
scale_from_unit <- function(scaled_vals, max_orig, min_orig) {
  vals_orig <- (scaled_vals * (max_orig - min_orig)) + min_orig
  return(vals_orig)
}
