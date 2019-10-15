#' Confidence Intervals for Stochastic Mediation Parameter Objects
#'
#' Compute confidence intervals for objects of class \code{medoutcon}, which
#' contain estimates produced by \code{medoutcon}.
#'
#' @param object An object of class \code{medoutcon}, as produced by invoking
#'  the function \code{tmle_medoutcon}, for which a confidence interval is to be
#'  computed.
#' @param parm A \code{numeric} vector indicating indices of \code{object$est}
#'  for which to return confidence intervals.
#' @param level A \code{numeric} indicating the level of the confidence interval
#'  to be computed.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats qnorm qlogis plogis
#' @importFrom assertthat assert_that
#'
#' @method confint medoutcon
#'
#' @export
#
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
    object$param %in% c("direct_effect", "indirect_effect")) {
    # NOTE: variance already scaled (i.e., Var(D)/n)
    se_eif <- sqrt(object$var)

    # compute the interval around the point estimate
    ci_theta <- ci_norm_bounds * se_eif + object$theta
  } else if (length(unique(object$outcome)) == 2 &&
    !(object$param %in% c("direct_effect", "indirect_effect"))) {
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

################################################################################

#' Summary for Stochastic Mediation Parameter Objects
#'
#' Print a convenient summary for objects of \code{S3} class \code{medoutcon}.
#'
#' @param object An object of class \code{medoutcon}, as produced by invoking
#'  the function \code{tmle_medoutcon}, for which a confidence interval is to be
#'  computed.
#' @param ... Other arguments. Not currently used.
#' @param ci_level A \code{numeric} indicating the level of the confidence
#'  interval to be computed.
#'
#' @importFrom stats confint
#'
#' @method summary medoutcon
#'
#' @export
#
summary.medoutcon <- function(object,
                              ...,
                              ci_level = 0.95) {
  # inference is currently limited to the efficient estimators
  if (object$type %in% c("onestep", "tmle")) {
    # compute confidence interval using the pre-defined method
    ci <- stats::confint(object, level = ci_level)

    # only print useful info about the mean of the efficient influence function
    eif_mean <- formatC(mean(object$eif), digits = 4, format = "e")

    # create output table from input object and confidence interval results
    out <- c(
      round(c(ci, object$var), digits = 4), eif_mean, object$type,
      object$param
    )
    names(out) <- c(
      "lwr_ci", "param_est", "upr_ci", "param_var", "eif_mean", "estimator",
      "param"
    )
  } else {
    out <- c(round(object$theta, digits = 6), object$type)
    names(out) <- c(
      "param_est", "estimator"
    )
  }
  print(noquote(out))
}

################################################################################

#' Print Method for Stochastic Mediation Parameter Objects
#'
#' The \code{print} method for objects of class \code{medoutcon}.
#'
#' @param x An object of class \code{medoutcon}.
#' @param ... Other options (not currently used).
#'
#' @method print medoutcon
#'
#' @export
#
print.medoutcon <- function(x, ...) {
  # inference is currently limited to the one-step efficient estimator
  # TODO: allow use for TML estimators once impelemented
  if (x$type %in% c("onestep", "tmle")) {
    print(x[c("theta", "var", "type", "param")])
  } else {
    print(x[c("theta", "type")])
  }
}

################################################################################

#' Bounding Numerical Precision
#'
#' Bounds extreme values to numerical (machine) precision, for use with
#' sensitive quantities like estimated propensity scores.
#'
#' @param vals A \code{numeric} vector of values in the interval [0, 1].
#'
#' @importFrom assertthat assert_that
#'
#' @keywords internal
#
bound_precision <- function(vals) {
  assertthat::assert_that(!(max(vals) > 1 | min(vals) < 0))
  vals[vals == 0] <- .Machine$double.neg.eps
  vals[vals == 1] <- 1 - .Machine$double.neg.eps
  return(vals)
}

################################################################################

#' Bounding Propensity Scores
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
#
bound_propensity <- function(vals, bounds = c(0.001, 0.999)) {
  assertthat::assert_that(!(max(vals) > 1 | min(vals) < 0))
  vals[vals < bounds[1]] <- bounds[1]
  vals[vals > bounds[2]] <- bounds[2]
  return(vals)
}

################################################################################

#' Scale Values to the Unit Interval [0, 1]
#'
#' @param vals A \code{numeric} vector of values to be scaled into the interval
#'  [0, 1].
#'
#' @keywords internal
#
scale_to_unit <- function(vals) {
  vals_scaled <- (vals - min(vals)) / (max(vals) - min(vals))
  return(vals_scaled)
}

################################################################################

#' Scale Values to the Original Scale
#'
#' @param scaled_vals A \code{numeric} vector of values scaled to lie in the
#'  interval [0, 1] by use of \code{\link{scale_to_unit}}.
#' @param max_orig The maximum of the values on the original scale.
#' @param min_orig The minimum of the values on the original scale.
#'
#' @keywords internal
#
scale_to_original <- function(scaled_vals, max_orig, min_orig) {
  vals_orig <- (scaled_vals * (max_orig - min_orig)) + min_orig
  return(vals_orig)
}
