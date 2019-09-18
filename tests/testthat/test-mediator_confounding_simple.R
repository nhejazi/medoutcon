context("Estimators show agreeable performance in simple example setting")

library(data.table)
library(stringr)
library(hal9001)
library(sl3)
library(caret)
library(SuperLearner)
set.seed(7128816)
source("eif_utils.R")
source("data_utils.R")

# 1) setup learners for the nuisance parameters
## caret hyperparameter-tuning model for random forest
SL.caretRF <- function(Y, X, newX, family, obsWeights, ...) {
  SL.caret(Y, X, newX, family, obsWeights, method = 'rf',  tuneLength = 3,
           trControl =  caret::trainControl(method = "LGOCV", search = 'random',
                                            verboseIter = TRUE), ...)
}
rf_caret_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.caretRF")

## caret hyperparameter-tuning model for xgboost
SL.caretXGB <- function(Y, X, newX, family, obsWeights, ...) {
  SL.caret(Y, X, newX, family, obsWeights, method = 'xgbTree', tuneLength = 3,
           trControl =  caret::trainControl(method = "LGOCV", search = 'random',
                                            verboseIter = TRUE), ...)
}
xgb_caret_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.caretXGB")

## instantiate learners and SL for continuous outcomes
mean_lrnr <- Lrnr_mean$new()
fglm_contin_lrnr <- Lrnr_glm_fast$new()
hal_contin_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5
)
stack_lrnrs_contin <- make_learner(Stack,
                                   mean_lrnr,
                                   fglm_contin_lrnr,
                                   hal_contin_lrnr,
                                   rf_caret_lrnr,
                                   xgb_caret_lrnr)
sl_contin_lrnr <- Lrnr_sl$new(learners = stack_lrnrs_contin,
                              metalearner = Lrnr_nnls$new(),
                              keep_extra = TRUE)

## instantiate learners and SL for binary outcomes
fglm_binary_lrnr <- Lrnr_glm_fast$new(family = binomial())
hal_binary_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5,
  family = "binomial"
)
logistic_metalearner <- make_learner(Lrnr_solnp, metalearner_logistic_binomial,
                                     loss_loglik_binomial)
stack_lrnrs_binary <- make_learner(Stack,
                                   mean_lrnr,
                                   fglm_binary_lrnr,
                                   hal_binary_lrnr,
                                   rf_caret_lrnr,
                                   xgb_caret_lrnr)
sl_binary_lrnr <- Lrnr_sl$new(learners = stack_lrnrs_binary,
                              metalearner = logistic_metalearner,
                              keep_extra = TRUE)

# 2) get data and column names for sl3 tasks (for convenience)
data <- sim_medoutcon_data()
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")

# 3) set up learners for nuisance parameters
## use SL including HAL for analyzing data
g_learners <- e_learners <- m_learners <- q_learners <- r_learners <-
  sl_binary_lrnr
u_learners <- v_learners <- sl_contin_lrnr
## use HAL by itself for testing functionality
#g_learners <- e_learners <- m_learners <- q_learners <- r_learners <-
  #hal_binary_lrnr
#u_learners <- v_learners <- hal_contin_lrnr

## test caret-based learners
#g_learners <- e_learners <- m_learners <- q_learners <- r_learners <-
  #xgb_caret_lrnr
#u_learners <- v_learners <- xgb_caret_lrnr
#g_learners <- e_learners <- m_learners <- q_learners <- r_learners <-
  #rf_caret_lrnr
#u_learners <- v_learners <- rf_caret_lrnr

# test different estimators
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
theta_os
