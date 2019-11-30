context("Estimators of parameter and EIF match manual analogs closely")
source("eif_utils.R")
source("data_utils.R")

# packages
library(data.table)
library(stringr)
library(tibble)
library(hal9001)
library(sl3)
# library(caret)
# library(SuperLearner)

# options
set.seed(7128816)
contrast <- c(0, 1)
aprime <- contrast[1]
astar <- contrast[2]
n_obs <- 1000

# 1) set up learners for nuisance parameters
if (FALSE) {
  ## caret hyperparameter-tuning model for random forest
  # SL.caretRF <- function(Y, X, newX, family, obsWeights, ...) {
  # SL.caret(Y, X, newX, family, obsWeights, method = 'rf',  tuneLength = 3,
  # trControl =  caret::trainControl(method = "LGOCV",
  # search = 'random',
  # verboseIter = TRUE), ...)
  # }
  # rf_caret_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.caretRF")

  ## caret hyperparameter-tuning model for xgboost
  # SL.caretXGB <- function(Y, X, newX, family, obsWeights, ...) {
  # SL.caret(Y, X, newX, family, obsWeights, method = 'xgbTree',
  # tuneLength = 3,
  # trControl =  caret::trainControl(method = "LGOCV",
  # search = 'random',
  # verboseIter = TRUE), ...)
  # }
  # xgb_caret_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.caretXGB")

  ## instantiate learners and SL for continuous outcomes
  mean_lrnr <- Lrnr_mean$new()
  fglm_contin_lrnr <- Lrnr_glm_fast$new()
  hal_contin_lrnr <- Lrnr_hal9001$new(
    fit_type = "glmnet", n_folds = 5
  )
  stack_lrnrs_contin <- make_learner(
    Stack,
    mean_lrnr,
    fglm_contin_lrnr,
    hal_contin_lrnr,
    rf_caret_lrnr,
    xgb_caret_lrnr
  )
  sl_contin_lrnr <- Lrnr_sl$new(
    learners = stack_lrnrs_contin,
    metalearner = Lrnr_nnls$new(),
    keep_extra = TRUE
  )

  ## instantiate learners and SL for binary outcomes
  fglm_binary_lrnr <- Lrnr_glm_fast$new(family = binomial())
  hal_binary_lrnr <- Lrnr_hal9001$new(
    fit_type = "glmnet", n_folds = 5,
    family = "binomial"
  )
  logistic_metalearner <- make_learner(
    Lrnr_solnp,
    metalearner_logistic_binomial,
    loss_loglik_binomial
  )
  stack_lrnrs_binary <- make_learner(
    Stack,
    mean_lrnr,
    fglm_binary_lrnr,
    hal_binary_lrnr,
    rf_caret_lrnr,
    xgb_caret_lrnr
  )
  sl_binary_lrnr <- Lrnr_sl$new(
    learners = stack_lrnrs_binary,
    metalearner = logistic_metalearner,
    keep_extra = TRUE
  )

  ## use SL including HAL for analyzing data
  g_learners <- e_learners <- m_learners <- q_learners <- r_learners <-
    sl_binary_lrnr
  u_learners <- v_learners <- sl_contin_lrnr
}

## use HAL by itself for testing functionality
hal_contin_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5
)
hal_binary_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5,
  family = "binomial"
)
g_learners <- e_learners <- m_learners <- q_learners <- r_learners <-
  hal_binary_lrnr
u_learners <- v_learners <- hal_contin_lrnr


# 3) get data and column names for sl3 tasks (for convenience)
data <- sim_medoutcon_data(n_obs = n_obs)
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")

# 4) test different estimators
theta_os <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = data$Z,
  M = data[, ..m_names], Y = data$Y,
  # effect = "direct",
  contrast = c(0, 1),
  g_learners = g_learners,
  e_learners = e_learners,
  m_learners = m_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  estimator = "onestep",
  estimator_args = list(cv_folds = 2)
)
summary(theta_os)

# 5) compute efficient influence function based on observed data
w <- as_tibble(data)[, w_names]
a <- data$A
z <- data$Z
m <- data$M
y <- data$Y

# compute parameter estimate and influence function with convenience functions
v <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) *
  pmaw(0, astar, w)
eif <- (a == aprime) / g(aprime, w) * pmaw(m, astar, w) /
  pm(m, z, aprime, w) * (y - my(m, z, aprime, w)) + (a == aprime) /
    g(aprime, w) * (u(z, w, aprime, astar) - intu(w, aprime, astar)) +
  (a == astar) / g(astar, w) * (intv(m, w, aprime) - v) + v
psi_os <- mean(eif)
var_eif <- var(eif) / n_obs

# 6) testing
test_that("Parameter estimate close to independent EIF estimates", {
  expect_equal(theta_os$theta, psi_os, tol = 0.04)
})

test_that("Variance estimate close to independent EIF variance", {
  expect_equal(theta_os$var, var_eif, tol = 0.01)
})

test_that("Mean of estimated EIF close to that of independent EIF", {
  expect_equal(abs(mean(theta_os$eif)), abs(mean(eif - psi_os)), tol = 0.001)
})
