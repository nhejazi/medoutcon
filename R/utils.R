#' Confidence Intervals for Stochastic Mediation Parameter Objects
#'
#' Compute confidence intervals for objects of class \code{medshift}, which
#' contain estimates produced by \code{medshift}.
#'
#' @param object An object of class \code{medshift}, as produced by invoking
#'  the function \code{tmle_medshift}, for which a confidence interval is to be
#'  computed.
#' @param parm A \code{numeric} vector indicating indices of \code{object$est}
#'  for which to return confidence intervals.
#' @param level A \code{numeric} indicating the level of the confidence interval
#'  to be computed.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats qnorm
#' @importFrom assertthat assert_that
#'
#' @method confint medshift
#'
#' @export
#
confint.medshift <- function(object,
                              parm = seq_len(object$psi),
                              level = 0.95,
                              ...) {
  # inference is currently limited to the one-step efficient estimator
  # TODO: allow use for TML estimators once impelemented
  assertthat::assert_that(object$type == "onestep")

  # first, let's get Z_(1 - alpha)
  ci_norm_bounds <- c(-1, 1) * abs(stats::qnorm(p = (1 - level) / 2))

  # compute the EIF variance multiplier for the CI
  # NOTE: the variance value is already scaled by length of observations
  se_eif <- sqrt(object$var)

  # compute the interval around the point estimate
  ci_theta <- ci_norm_bounds * se_eif + object$theta

  # set up output CI object
  ci_out <- c(ci_theta[1], object$theta, ci_theta[2])
  names(ci_out) <- c("lwr_ci", "est", "upr_ci")
  return(ci_out)
}

################################################################################

#' Summary for Stochastic Mediation Parameter Objects
#'
#' Print a convenient summary for objects of \code{S3} class \code{medshift}.
#'
#' @param object An object of class \code{medshift}, as produced by invoking
#'  the function \code{tmle_medshift}, for which a confidence interval is to be
#'  computed.
#' @param ... Other arguments. Not currently used.
#' @param ci_level A \code{numeric} indicating the level of the confidence
#'  interval to be computed.
#'
#' @importFrom stats confint
#'
#' @method summary medshift
#'
#' @export
#
summary.medshift <- function(object,
                              ...,
                              ci_level = 0.95) {
  # inference is currently limited to the one-step efficient estimator
  if (object$type %in% c("onestep", "tmle")) {
    # compute confidence interval using the pre-defined method
    ci <- stats::confint(object, level = ci_level)

    # only print useful info about the mean of the efficient influence function
    eif_mean <- formatC(mean(object$eif), digits = 4, format = "e")

    # create output table from input object and confidence interval results
    out <- c(round(c(ci, object$var), digits = 4), eif_mean, object$type)
    names(out) <- c(
      "lwr_ci", "param_est", "upr_ci", "param_var", "eif_mean", "estimator"
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
#' The \code{print} method for objects of class \code{medshift}.
#'
#' @param x An object of class \code{medshift}.
#' @param ... Other options (not currently used).
#'
#' @method print medshift
#'
#' @export
#
print.medshift <- function(x, ...) {
  # inference is currently limited to the one-step efficient estimator
  # TODO: allow use for TML estimators once impelemented
  if (x$type == "one-step efficient") {
    print(x[c("theta", "var", "type")])
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
  assertthat::assert_that(!(max(vals) >= 1 | min(vals) <= 0))
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
bound_propensity <- function(vals, bounds = c(0.01, 0.99)) {
  assertthat::assert_that(!(max(vals) >= 1 | min(vals) <= 0))
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

################################################################################

#' Numerical Integration for Weighted Treatment Mechanism
#'
#' In the case of modified treatment policies, it is necessary to numerically
#' evaluate an integral over the domain of the treatment mechanism. This is a
#' simple procedure to numerically compute such an integral based on Monte
#' Carlo importance sampling from a uniform distribution.
#'
#' @param g_mech The estimated conditional density corresponding to the natural
#'  or shifted value of the treatment for modified treatment policies based on
#'  the observed values of the treatment.
#' @param a_vals The observed values of the treatment used in computing the
#'  conditional density estimates provided in the argument \code{g_mech}.
#' @param weighting A \code{numeric} vector of weights corresponding to various
#'  regression functions estiated based on the observed values of the treatment.
#'  Such values correspond to a multiplicative weight applied to the estimated
#'  conditional density.
#' @param int_grid_prop A \code{numeric} scalar corresponding to the propotion
#'  of points to be used in numerically evaluating an integral over the domain
#'  of the natural/shifted conditional density of the treatment. This is only
#'  relevant in the case of modified treatment policies for continuous-valued
#'  exposures; it is irrelevant for incremental propensity score interventions.
#'
#' @keywords internal
#
integrate_over_g <- function(g_mech, a_vals, weighting, int_grid_prop = 0.5) {
  # numerical integration over the domain of A via Monte Carlo
  min_a_shifted <- min(a_vals)
  max_a_shifted <- max(a_vals)
  range_a_shifted <- max_a_shifted - min_a_shifted

  # sample points uniformly distributed over the treatment mechanism
  int_grid_size <- round(length(a_vals) * int_grid_prop)
  int_grid_points <- sample(length(a_vals), int_grid_size)
  g_mech_grid <- g_mech[int_grid_points]

  # choose same grid of uniformly sampled points for the weighting component
  weighting_grid <- weighting[int_grid_points]

  # compute weighted combination of terms for numerical integration
  integ_mc <- (range_a_shifted / int_grid_size) * g_mech_grid * weighting_grid
  return(integ_mc)
}
