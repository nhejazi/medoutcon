context("Estimators of nuisance parameters match manual analogs closely")
source("utils_interventional.R")

# packages
library(data.table)
library(stringr)
library(tibble)
library(hal9001)
library(sl3)

# options
contrast <- c(0, 1)
aprime <- contrast[1]
astar <- contrast[2]
set.seed(27158)
n_samp <- 5000

# set up learners for each nuisance parameter
hal_binomial_lrnr <- Lrnr_hal9001$new(
  family = "binomial",
  fit_control = list(
    n_folds = 5,
    use_min = TRUE,
    type.measure = "mse",
    lambda.min.ratio = 1 / n_samp
  )
)
hal_gaussian_lrnr <- Lrnr_hal9001$new(
  family = "gaussian",
  fit_control = list(
    n_folds = 5,
    use_min = TRUE,
    type.measure = "mse",
    lambda.min.ratio = 1 / n_samp
  )
)
g_learners <- h_learners <- b_learners <- q_learners <- r_learners <- hal_binomial_lrnr
u_learners <- v_learners <- hal_gaussian_lrnr

# simulate smaller data set for computing estimates
data <- sim_medoutcon_data(n_obs = n_samp)
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")
data[, obs_weights := 1]
w <- as_tibble(data)[, w_names]
a <- data$A
z <- data$Z
m <- data$M
y <- data$Y

# compute estimates of nuisance parameters
## fit propensity score
g_out <- fit_treat_mech(
  train_data = data, valid_data = data, contrast = c(aprime, astar),
  learners = g_learners, w_names = w_names, type = "g"
)
test_that("Estimates of propensity score are close to the truth", {
  expect_equal(g_out$treat_est_train$treat_pred_A_star, g(astar, w),
    tol = 0.05
  )
})
test_that("MSE of propensity score estimates is sufficiently low", {
  g_mse <- mean((g_out$treat_est_train$treat_pred_A_star - g(astar, w))^2)
  expect_lt(g_mse, 1e-3)
})

## fit propensity score conditioning on mediators
h_out <- fit_treat_mech(
  train_data = data, valid_data = data, contrast = c(aprime, astar),
  learners = h_learners, w_names = w_names,
  m_names = m_names, type = "h"
)
test_that("Estimates of mediator propensity score are close to the truth", {
  expect_equal(h_out$treat_est_train$treat_pred_A_prime, e(aprime, m, w),
    tol = 0.025
  )
})
test_that("MSE of mediator propensity score estimates is sufficiently low", {
  h_mse <- mean((h_out$treat_est_train$treat_pred_A_prime - e(aprime, m, w))^2)
  expect_lt(h_mse, 1e-3)
})

## fit outcome regression
b_out <- fit_out_mech(
  train_data = data, valid_data = data, contrast = c(aprime, astar),
  learners = b_learners, m_names = m_names, w_names = w_names
)
test_that("Estimates of outcome regression are close to the truth", {
  expect_equal(b_out$b_est_train$b_pred_A_prime, my(m, z, aprime, w),
    tol = 0.15
  )
})
test_that("MSE of outcome regression estimates is sufficiently low", {
  b_mse <- mean((b_out$b_est_train$b_pred_A_prime - my(m, z, aprime, w))^2)
  expect_lt(b_mse, 0.01)
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
test_that("Estimates of confounder regression q are close to the truth", {
  expect_equal(q_out$moc_est_train_Z_one$moc_pred_A_prime, pz(1, aprime, w),
    tol = 0.05
  )
})
test_that("MSE of confounder regression q estimates is sufficiently low", {
  q_mse <- mean((q_out$moc_est_train_Z_one$moc_pred_A_prime -
    pz(1, aprime, w))^2)
  expect_lt(q_mse, 1e-3)
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
test_that("Estimates of confounder regression r are close to the truth", {
  expect_equal(r_out$moc_est_train_Z_one$moc_pred_A_prime,
    r(1, aprime, m, w),
    tol = 0.075
  )
})
test_that("MSE of confounder regression r estimates is sufficiently low", {
  r_mse <- mean((r_out$moc_est_train_Z_one$moc_pred_A_prime -
    r(1, aprime, m, w))^2)
  expect_lt(r_mse, 2e-3)
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
test_that("Estimates of pseudo-outcome regression are close to the truth", {
  expect_equal(u_out$u_pred, u(z, w, aprime, astar), tol = 0.12)
})
test_that("MSE of pseudo-outcome regression estimates is sufficiently low", {
  u_mse <- mean((u_out$u_pred - u(z, w, aprime, astar))^2)
  expect_lt(u_mse, 0.01)
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
v_star <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) *
  pmaw(0, astar, w)

test_that("Estimates of pseudo-outcome regression are close to the truth", {
  expect_equal(v_out$v_pred, v_star, tol = 0.08)
})
test_that("MSE of pseudo-outcome regression estimates is sufficiently low", {
  v_mse <- mean((v_out$v_pred - v_star)^2)
  expect_lt(v_mse, 0.005)
})
test_that("Estimates of pseudo-outcome used in v are close to the truth", {
  expect_equal(v_out$v_pseudo, intv(m, w, aprime), tol = 0.1)
})
test_that("MSE of pseudo-outcome used in v estimation is sufficiently low", {
  v_pseudo_mse <- mean((v_out$v_pseudo - intv(m, w, aprime))^2)
  expect_lt(v_pseudo_mse, 0.005)
})
