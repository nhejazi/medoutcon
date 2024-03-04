################################################################################
# Unit Tests
################################################################################
context("Two-phase sampling EIF convenience function")

test_that("two_phase_eif returns an uncentered EIF", {
  # generate fake inputs
  R <- c(1, 0, 0, 0, 1, 1)
  two_phase_weights <- c(rep(1 / 2, 3), rep(2, 3))
  partial_eif <- c(1 / 2, rep(-1 / 2, 2))
  eif_predictions <- c(rep(1 / 3, 3), rep(-1 / 3, 3))
  plugin_est <- 1

  # independently compute the adjusted EIF
  uncentered_eif <- R * two_phase_weights *
    c(partial_eif[1], 0, 0, 0, partial_eif[2], partial_eif[3]) +
    (1 - R * two_phase_weights) * eif_predictions + plugin_est

  # assert that result is identical to two_phase_eif's
  expect_equal(
    two_phase_eif(
      R = R,
      two_phase_weights = two_phase_weights,
      eif = partial_eif,
      eif_predictions = eif_predictions,
      plugin_est = plugin_est
    ),
    uncentered_eif
  )
})


context(paste(
  "Estimators of nuisance parameters match manual analogs closely",
  "in two-phase sampling designs"
))
source("utils_interventional.R")
source("utils_natural.R")

# packages
library(data.table)
library(stringr)
library(tibble)
library(hal9001)
library(sl3)
library(origami)

# options
contrast <- c(0, 1)
aprime <- contrast[1]
astar <- contrast[2]
set.seed(245627)
n_samp <- 5000

# set up learners for each nuisance parameter
hal_binomial_lrnr <- Lrnr_hal9001$new(
  family = "binomial",
  fit_control = list(
    n_folds = 5,
    use_min = TRUE
  )
)
hal_gaussian_lrnr <- Lrnr_hal9001$new(
  family = "gaussian",
  fit_control = list(
    n_folds = 5,
    use_min = TRUE
  )
)
g_learners <- h_learners <- b_learners <- q_learners <- r_learners <-
  hal_binomial_lrnr
u_learners <- v_learners <- hal_gaussian_lrnr

# simulate smaller data set for computing estimates
data <- sim_medoutcon_data(n_obs = n_samp)
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")
data[, `:=`(
  R = rbinom(n_samp, 1, 0.9),
  two_phase_weights = 1,
  obs_weights = 1
)]
w <- as_tibble(data)[, w_names]
a <- data$A
z <- data$Z
m <- data$M
y <- data$Y
R <- data$R

# compute estimates of nuisance parameters
## fit propensity score
g_out <- fit_treat_mech(
  train_data = data, valid_data = data, contrast = c(aprime, astar),
  learners = g_learners, w_names = w_names, type = "g"
)
test_that("MSE of propensity score estimates is sufficiently low", {
  g_mse <- mean((g_out$treat_est_train$treat_pred_A_star - g(astar, w))^2)
  expect_lt(g_mse, 0.05)
})

## fit propensity score conditioning on mediators
h_out <- fit_treat_mech(
  train_data = data, valid_data = data, contrast = c(aprime, astar),
  learners = h_learners, w_names = w_names,
  m_names = m_names, type = "h"
)
test_that("MSE of mediator propensity score estimates is sufficiently low", {
  h_mse <- mean(
    (
      h_out$treat_est_train$treat_pred_A_prime -
        e(aprime, m[data$R == 1], w[data$R == 1, ])
    )^2
  )
  expect_lt(h_mse, 0.01)
})

## fit outcome regression
b_out <- fit_out_mech(
  train_data = data, valid_data = data, contrast = c(aprime, astar),
  learners = b_learners, m_names = m_names, w_names = w_names
)
test_that("MSE of outcome regression estimates is sufficiently low", {
  b_mse <- mean(
    (b_out$b_est_train$b_pred_A_prime -
      my(m[data$R == 1], z[data$R == 1], aprime, w[data$R == 1, ])
    )^2
  )
  expect_lt(b_mse, 0.02)
})

## fit mediator-outcome confounder regression, excluding mediator(s)
q_out <- fit_moc_mech(
  train_data = data,
  valid_data = data,
  contrast = contrast,
  learners = q_learners,
  m_names = m_names,
  w_names = w_names,
  type = "q"
)
test_that("MSE of confounder regression q estimates is sufficiently low", {
  q_mse <- mean((q_out$moc_est_train_Z_one$moc_pred_A_prime -
    pz(1, aprime, w))^2)
  expect_lt(q_mse, 0.01)
})

## fit mediator-outcome confounder regression, conditioning on mediator(s)
r_out <- fit_moc_mech(
  train_data = data,
  valid_data = data,
  contrast = contrast,
  learners = r_learners,
  m_names = m_names,
  w_names = w_names,
  type = "r"
)
test_that("MSE of confounder regression r estimates is sufficiently low", {
  r_mse <- mean(
    (
      r_out$moc_est_train_Z_one$moc_pred_A_prime -
        r(1, aprime, m[data$R == 1], w[data$R == 1, ])
    )^2
  )
  expect_lt(r_mse, 0.01)
})

# data for fitting nuisance parameters with pseudo-outcomes
data_a_prime <- data.table::copy(data)[, A := aprime]
data_a_star <- data.table::copy(data)[, A := astar]

# fit u
u_out <- fit_nuisance_u(
  train_data = data,
  valid_data = data_a_prime,
  learners = u_learners,
  b_out = b_out,
  q_out = q_out,
  r_out = r_out,
  g_out = g_out,
  h_out = h_out,
  w_names = w_names
)
test_that("MSE of pseudo-outcome regression estimates is sufficiently low", {
  u_mse <- mean(
    (u_out$u_pred -
      u(z[data$R == 1], w[data$R == 1, ], aprime, astar)
    )^2
  )
  expect_lt(u_mse, 0.02)
})

## fit v
v_out <- fit_nuisance_v(
  train_data = data,
  valid_data = data_a_star,
  contrast = contrast,
  learners = v_learners,
  b_out = b_out,
  q_out = q_out,
  m_names = m_names,
  w_names = w_names
)
v_star <- intv(1, w[data$R == 1, ], aprime) * pmaw(1, astar, w[data$R == 1, ]) +
  intv(0, w[data$R == 1, ], aprime) * pmaw(0, astar, w[data$R == 1, ])

test_that("MSE of pseudo-outcome regression estimates is sufficiently low", {
  v_mse <- mean((v_out$v_pred - v_star)^2)
  expect_lt(v_mse, 0.01)
})
test_that("MSE of pseudo-outcome used in v estimation is sufficiently low", {
  v_pseudo_mse <- mean(
    (v_out$v_pseudo - intv(m[data$R == 1], w[data$R == 1, ], aprime))^2
  )
  expect_lt(v_pseudo_mse, 0.01)
})

################################################################################
# Integration tests for estimators
################################################################################

context(paste(
  "Estimators of natural effects are close to true parameters in",
  "two-phase sampling designs"
))
source("utils_natural.R")

# packages
library(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(hal9001)
library(glmnet)
library(sl3)
library(SuperLearner)

# options
set.seed(8123421)
n_obs <- 500

# 1) get data and column names for sl3 tasks (for convenience)
data <- make_nide_data(n_obs = n_obs)
R <- rbinom(n_obs, 1, 0.9)
R[data$Y > 5] <- 1
two_phase_weights <- rep(1 / 0.9, nrow(data))
two_phase_weights[data$Y > 5] <- 1
data[, `:=`(
  R = R,
  two_phase_weights = two_phase_weights,
  obs_weights = 1
)]
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "Z")

# 2) use simpler SLs for testing functionality
mean_lrnr <- Lrnr_mean$new()
fglm_lrnr <- Lrnr_glm_fast$new()
lasso_lrnr <- Lrnr_glmnet$new(alpha = 1, nfolds = 3L)
enet_lrnr <- Lrnr_glmnet$new(alpha = 0.5, nfolds = 3L)
rf_lrnr <- Lrnr_ranger$new(
  num.trees = 1000, sample.fraction = 0.7,
  oob.error = FALSE
)
logistic_meta <- Lrnr_solnp$new(
  metalearner_logistic_binomial,
  loss_loglik_binomial
)
sl_binary <- Lrnr_sl$new(
  learners = list(
    rf_lrnr,
    fglm_lrnr,
    lasso_lrnr,
    enet_lrnr,
    mean_lrnr
  ),
  metalearner = logistic_meta
)
sl_contin <- Lrnr_sl$new(
  learners = list(
    rf_lrnr,
    fglm_lrnr,
    lasso_lrnr,
    enet_lrnr,
    mean_lrnr
  ),
  metalearner = Lrnr_nnls$new()
)

## nuisance functions with data components have binary outcomes
g_learners <- h_learners <- q_learners <- r_learners <- rf_lrnr

## nuisance functions with pseudo-outcomes have continuous outcomes
d_learners <- u_learners <- v_learners <- b_learners <- rf_lrnr


# 3) test different estimators
nde_os <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = NULL,
  M = data[, ..m_names], Y = data$Y, R = data$R,
  two_phase_weights = data$two_phase_weights,
  g_learners = g_learners,
  h_learners = h_learners,
  b_learners = b_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  d_learners = d_learners,
  effect = "direct",
  estimator = "onestep",
  estimator_args = list(cv_folds = 5, max_iter = 0, tiltmod_tol = 5)
)
summary(nde_os)

nie_os <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = NULL,
  M = data[, ..m_names], Y = data$Y, R = data$R,
  two_phase_weights = data$two_phase_weights,
  g_learners = g_learners,
  h_learners = h_learners,
  b_learners = b_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  d_learners = d_learners,
  effect = "indirect",
  estimator = "onestep",
  estimator_args = list(cv_folds = 5, max_iter = 0, tiltmod_tol = 5)
)
summary(nie_os)

nde_tmle <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = NULL,
  M = data[, ..m_names], Y = data$Y, R = data$R,
  two_phase_weights = data$two_phase_weights,
  g_learners = g_learners,
  h_learners = h_learners,
  b_learners = b_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  d_learners = d_learners,
  effect = "direct",
  estimator = "tmle",
  estimator_args = list(cv_folds = 5, max_iter = 10, tiltmod_tol = 5)
)
summary(nde_tmle)

nie_tmle <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = NULL,
  M = data[, ..m_names], Y = data$Y, R = data$R,
  two_phase_weights = data$two_phase_weights,
  g_learners = g_learners,
  h_learners = h_learners,
  b_learners = b_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  d_learners = d_learners,
  effect = "indirect",
  estimator = "tmle",
  estimator_args = list(cv_folds = 5, max_iter = 10, tiltmod_tol = 5)
)
summary(nie_tmle)


# 4) compute approximate truth and efficiency bound
sim_truth <- get_truth_nide(
  n_obs = 1e6, binary_outcome = FALSE, EIC = TRUE
)
EY_A1_Z1 <- sim_truth$EY_A1_Z1
EY_A1_Z0 <- sim_truth$EY_A1_Z0
EY_A0_Z1 <- sim_truth$EY_A0_Z1
EY_A0_Z0 <- sim_truth$EY_A0_Z0
nde_true <- mean(EY_A1_Z0 - EY_A0_Z0)
nie_true <- mean(EY_A1_Z1 - EY_A1_Z0)


# 5) testing estimators for the NDE
test_that("NDE: One-step estimate is near DGP truth", {
  expect_equal(nde_os$theta, nde_true,
    tol = 1.96 * sqrt(var(nde_os$eif) / n_obs)
  )
})

test_that("NDE: TML estimate is near DGP truth", {
  expect_equal(nde_tmle$theta, nde_true,
    tol = 1.96 * sqrt(var(nde_tmle$eif) / n_obs)
  )
})

test_that("NDE: Mean of estimated EIF is nearly zero for the one-step", {
  expect_lt(abs(mean(nde_os$eif)), 1e-15)
})

# 6) testing estimators for the NIE
test_that("NIE: One-step estimate is near DGP truth", {
  expect_equal(nie_os$theta, nie_true,
    tol = 1.96 * sqrt(var(nie_os$eif) / n_obs)
  )
})

test_that("NIE: TML estimate is near DGP truth", {
  expect_equal(nie_tmle$theta, nie_true,
    tol = 1.96 * sqrt(var(nie_tmle$eif) / n_obs)
  )
})

test_that("NIE: Mean of estimated EIF is nearly zero for the one-step", {
  expect_lt(abs(mean(nie_os$eif)), 1e-15)
})
