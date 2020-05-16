#' Custom Highly Adaptive Lasso for SuperLearner
#'
#' Wrapper for \pkg{SuperLearner} for a customized Highly Adaptive Lasso
#'
#' @param Y A \code{numeric} of outcomes.
#' @param X A \code{matrix} of predictors/covariates.
#' @param newX A matrix of new observations on which to obtain predictions. The
#'  default of \code{NULL} computes predictions on training inputs \code{X}.
#' @param family Not used by the function directly, but meant to ensure
#'  compatibility with \code{SuperLearner}.
#' @param obsWeights Not used by the function directly, but meant to ensure
#'  compatibility with \code{SuperLearner}. These are passed to
#'  \code{\link[glmnet]{cv.glmnet}}.
#' @param id A \code{numeric} of observation IDs.
#' @param alpha A \code{numeric} indicating the mixing parameter for the Lasso
#'  (L1) and Ridge (L2) penalties. The default of \code{alpha = 1} is for the
#'  Lasso L1 penalty.
#' @param nfolds Integer for the number of folds to be used when splitting the
#'  data for cross-validation. This defaults to 5.
#' @param nlambda A \code{numeric} corresponding to the length of the sequence
#'  to be used for selecting an optimal value of the penalization parameter.
#' @param useMin Determines which lambda is selected from \code{cv.glmnet}.
#'  \code{TRUE} corresponds to \code{"lambda.min"} and \code{FALSE} corresponds
#'  to \code{"lambda.1se"}.
#' @param loss A \code{character} indicating the loss function that ought to be
#'  considered. Choices \code{"mse"} for mean squared error, \code{"deviance"},
#'  and any others available through \code{\link[glmnet]{cv.glmnet}}.
#' @param ... Placeholder (ignored).
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats predict as.formula model.matrix
#'
#' @export
#'
#' @return An object of class \code{SL.halglmnet} with a fitted \code{glmnet}
#'  object and corresponding predictions based on the input data.
SL.halglmnet <- function(Y, X, newX, family, obsWeights,
                         id, alpha = 1, nfolds = 5, nlambda = 1000,
                         useMin = TRUE, loss = "mse", ...) {
  # NOTE: `family` is ignored to allow use of `lambda.min.ratio`
  if (!is.matrix(X)) {
    formula <- stats::as.formula(paste("~ -1 + .^", eval(ncol(X))))
    X <- stats::model.matrix(formula, X)
    newX <- stats::model.matrix(formula, newX)
  }
  fitCV <- glmnet::cv.glmnet(
    x = X, y = Y, weights = obsWeights,
    lambda = NULL, type.measure = loss,
    alpha = alpha, nfolds = nfolds,
    family = "gaussian",
    nlambda = nlambda,
    lambda.min.ratio = 1 / nrow(X),
    ...
  )
  pred <- stats::predict(fitCV,
    newx = newX, type = "response",
    s = ifelse(useMin, "lambda.min", "lambda.1se")
  )
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.halglmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

################################################################################

#' predict.SL.halglmnet
#'
#' Predict method for objects of class \code{SL.halglmnet}
#'
#' @param object A fitted object of class \code{glmnet}.
#' @param newdata A matrix of new observations on which to obtain predictions.
#' @param remove_extra_cols A \code{logical} indicating whether extra columns
#'  should be removed.
#' @param add_missing_cols A \code{logical} indicating whether missing columns
#'  should be added to make the \code{newdata} and training data consistent.
#' @param ... Placeholder (ignored).
#'
#' @importFrom stats predict as.formula model.matrix
#'
#' @export
#'
#' @return A \code{numeric} vector of predictions from a \code{SL.halglmnet}
#'  object based on the provided \code{newdata}.
predict.SL.halglmnet <- function(object, newdata, remove_extra_cols = FALSE,
                                 add_missing_cols = FALSE, ...) {
  if (!is.matrix(newdata)) {
    formula <- stats::as.formula(paste("~ -1 + .^", eval(ncol(newdata))))
    newdata <- stats::model.matrix(formula, newdata)
  }
  original_cols <- rownames(object$object$glmnet.fit$beta)
  if (remove_extra_cols) {
    extra_cols <- setdiff(colnames(newdata), original_cols)
    if (length(extra_cols) > 0) {
      warning(paste(
        "Removing extra columns in prediction data:",
        paste(extra_cols, collapse = ", ")
      ))
      newdata <- newdata[, !colnames(newdata) %in% extra_cols,
        drop = FALSE
      ]
    }
  }
  if (add_missing_cols) {
    missing_cols <- setdiff(original_cols, colnames(newdata))
    if (length(missing_cols) > 0) {
      warning(paste(
        "Adding missing columns in prediction data:",
        paste(missing_cols, collapse = ", ")
      ))
      new_cols <- matrix(0,
        nrow = nrow(newdata),
        ncol = length(missing_cols)
      )
      colnames(new_cols) <- missing_cols
      newdata <- cbind(newdata, new_cols)
      newdata <- newdata[, original_cols]
    }
  }
  pred <- stats::predict(object$object,
    newx = newdata, type = "response",
    s = ifelse(object$useMin, "lambda.min",
      "lambda.1se"
    )
  )
  return(pred)
}
