context("Estimators of parameter and EIF match manual analogs closely")
source("utils_eif.R")
source("utils_data.R")

# packages
library(data.table)
library(stringr)
library(tibble)
library(hal9001)
library(sl3)
library(SuperLearner)

# options
set.seed(27158)
contrast <- c(0, 1)
aprime <- contrast[1]
astar <- contrast[2]
n_obs <- 5000

# 1) get data and column names for sl3 tasks (for convenience)
data <- sim_medoutcon_data(n_obs = n_obs)
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")

# 2) use custom HAL and SL for testing functionality
bounding_lrnr <- Lrnr_bound$new(bound = 1e-6)
mean_lrnr <- Lrnr_mean$new()
hal_custom_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.halglmnet")
hal_gaussian_lrnr <- Lrnr_hal9001$new(family = "gaussian",
                                      type.measure = "mse", n_folds = 5,
                                      use_min = FALSE, max_degree = NULL,
                                      lambda.min.ratio = 1 / n_obs)
hal_bounded_lrnr <- Pipeline$new(hal_gaussian_lrnr, bounding_lrnr)
sl <- Lrnr_sl$new(learners = list(hal_custom_lrnr,
                                  hal_bounded_lrnr,
                                  mean_lrnr),
                  metalearner = Lrnr_nnls$new())

## nuisance functions with data components as outcomes
g_learners <- e_learners <- m_learners <- q_learners <- r_learners <- sl

## nuisance functions with pseudo-outcomes need Gaussian HAL
u_learners <- v_learners <- hal_gaussian_lrnr

# 3) test different estimators
theta_os <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = data$Z,
  M = data[, ..m_names], Y = data$Y,
  contrast = c(0, 1),
  g_learners = g_learners,
  e_learners = e_learners,
  m_learners = m_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  estimator = "onestep",
  estimator_args = list(cv_folds = 5, max_iter = 0, tiltmod_tol = 5)
)
summary(theta_os)

theta_tmle <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = data$Z,
  M = data[, ..m_names], Y = data$Y,
  contrast = c(0, 1),
  g_learners = g_learners,
  e_learners = e_learners,
  m_learners = m_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  estimator = "tmle",
  estimator_args = list(cv_folds = 5, max_iter = 5, tiltmod_tol = 5)
)
summary(theta_tmle)

# 4) compute estimate and influence function with convenience functions
w <- as_tibble(data)[, w_names]
a <- data$A
z <- data$Z
m <- data$M
y <- data$Y
v <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) *
  pmaw(0, astar, w)
eif <- (a == aprime) / g(aprime, w) * pmaw(m, astar, w) /
  pm(m, z, aprime, w) * (y - my(m, z, aprime, w)) + (a == aprime) /
    g(aprime, w) * (u(z, w, aprime, astar) - intu(w, aprime, astar)) +
  (a == astar) / g(astar, w) * (intv(m, w, aprime) - v) + v
psi_indep <- mean(v)
var_indep <- var(eif) / n_obs

# 5) testing one-step estimator
test_that("One-step estimate close to independent EIF estimates", {
  expect_equal(theta_os$theta, psi_indep, tol = 0.005)
})

test_that("EIF variance of one-step is close to independent EIF variance", {
  expect_equal(theta_os$var, var_indep, tol = 0.0001)
})

test_that("Mean of estimated EIF is nearly zero for the one-step", {
  expect_lt(abs(mean(theta_os$eif)), 1e-7)
})

# 6) testing TML estimator
test_that("TML estimate close to independent EIF estimates", {
  expect_equal(theta_tmle$theta, psi_indep, tol = 0.005)
})

test_that("EIF variance of TMLE is close to independent EIF variance", {
  expect_equal(theta_tmle$var, var_indep, tol = 0.0001)
})

test_that("Mean of estimated EIF is nearly zero for the TMLE", {
  expect_lt(abs(mean(theta_tmle$eif)), 1e-7)
})
