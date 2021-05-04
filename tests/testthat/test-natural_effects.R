context("Estimators of natural effects match manual analogs closely")
source("utils_natural.R")

# packages
library(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(hal9001)
library(sl3)
library(SuperLearner)

# options
set.seed(27158)
n_obs <- 2000

# 1) get data and column names for sl3 tasks (for convenience)
data <- make_nide_data(n_obs = n_obs)
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "Z")

# 2) use custom HAL and SL for testing functionality
bounding_lrnr <- Lrnr_bound$new(bound = 1e-6)
fglm_lrnr <- Lrnr_glm_fast$new(family = binomial())
mean_lrnr <- Lrnr_mean$new()
hal_gaussian_lrnr <- Lrnr_hal9001$new(
  family = "gaussian",
  fit_control = list(
    n_folds = 5,
    use_min = TRUE,
    type.measure = "mse",
    lambda.min.ratio = 1 / n_obs
  )
)
hal_binomial_lrnr <- Lrnr_hal9001$new(
  family = "binomial",
  fit_control = list(
    n_folds = 5,
    use_min = TRUE,
    type.measure = "mse",
    lambda.min.ratio = 1 / n_obs
  )
)
hal_bounded_lrnr <- Pipeline$new(hal_gaussian_lrnr, bounding_lrnr)
sl <- Lrnr_sl$new(
  learners = list(
    hal_bounded_lrnr,
    hal_binomial_lrnr,
    fglm_lrnr,
    mean_lrnr
  ),
  metalearner = Lrnr_nnls$new()
)

## nuisance functions with data components as outcomes
g_learners <- h_learners <- q_learners <- r_learners <- sl

## nuisance functions with pseudo-outcomes need Gaussian HAL
u_learners <- v_learners <- b_learners <- hal_gaussian_lrnr

# 3) test different estimators
nde_os <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = NULL,
  M = data[, ..m_names], Y = data$Y,
  g_learners = g_learners,
  h_learners = h_learners,
  b_learners = b_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  effect = "direct",
  estimator = "onestep",
  estimator_args = list(cv_folds = 5, max_iter = 0, tiltmod_tol = 5)
)
summary(nde_os)

nie_os <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = NULL,
  M = data[, ..m_names], Y = data$Y,
  g_learners = g_learners,
  h_learners = h_learners,
  b_learners = b_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  effect = "indirect",
  estimator = "onestep",
  estimator_args = list(cv_folds = 5, max_iter = 0, tiltmod_tol = 5)
)
summary(nie_os)

nde_tmle <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = NULL,
  M = data[, ..m_names], Y = data$Y,
  g_learners = g_learners,
  h_learners = h_learners,
  b_learners = b_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  effect = "direct",
  estimator = "tmle",
  estimator_args = list(cv_folds = 5, max_iter = 5, tiltmod_tol = 5)
)
summary(nde_tmle)

nie_tmle <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = NULL,
  M = data[, ..m_names], Y = data$Y,
  g_learners = g_learners,
  h_learners = h_learners,
  b_learners = b_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  effect = "indirect",
  estimator = "tmle",
  estimator_args = list(cv_folds = 5, max_iter = 5, tiltmod_tol = 5)
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
var_nde_eff <- var(sim_truth$EIC_NDE) / n_obs
var_nie_eff <- var(sim_truth$EIC_NIE) / n_obs


# 5) testing estimators for the NDE
test_that("NDE: One-step estimate is near approximate truth", {
  expect_equal(nde_os$theta, nde_true, tol = 0.05)
})

test_that("NDE: TML estimate is near approximate truth", {
  expect_equal(nde_tmle$theta, nde_true, tol = 0.05)
})

test_that("NDE: EIF variance of one-step is near the true EIF variance", {
  expect_equal(nde_os$var, var_nde_eff, tol = 0.02)
})

test_that("NDE: EIF variance of TMLE is near the true EIF variance", {
  expect_equal(nde_tmle$var, var_nde_eff, tol = 0.02)
})

test_that("NDE: Mean of estimated EIF is nearly zero for the one-step", {
  expect_lt(abs(mean(nde_os$eif)), 1e-15)
})

test_that("NDE: Mean of estimated EIF is nearly zero for the TMLE", {
  expect_lt(abs(mean(nde_tmle$eif)), 1e-10)
})


# 6) testing estimators for the NIE
test_that("NIE: One-step estimate is near approximate truth", {
  expect_equal(nie_os$theta, nie_true, tol = 0.05)
})

test_that("NIE: TML estimate is near approximate truth", {
  expect_equal(nie_tmle$theta, nie_true, tol = 0.05)
})

test_that("NIE: EIF variance of one-step is near the true EIF variance", {
  expect_equal(nie_os$var, var_nie_eff, tol = 0.03)
})

test_that("NIE: EIF variance of TMLE is near the true EIF variance", {
  expect_equal(nie_tmle$var, var_nie_eff, tol = 0.03)
})

test_that("NIE: Mean of estimated EIF is nearly zero for the one-step", {
  expect_lt(abs(mean(nie_os$eif)), 1e-15)
})

test_that("NIE: Mean of estimated EIF is nearly zero for the TMLE", {
  expect_lt(abs(mean(nie_tmle$eif)), 1e-10)
})
