context("Estimators of nuisance parameters match EIF-based analogs")

# packages and options
library(data.table)
library(stringr)
library(tibble)
library(hal9001)
library(sl3)
source("eif_utils.R")
source("data_utils.R")
contrast <- c(0, 1)
aprime <- contrast[1]
astar <- contrast[2]
set.seed(7128816)


# set up learners for each nuisance parameter
g_learners <- e_learners <- m_learners <- q_learners <- r_learners <-
  Lrnr_hal9001$new(family = "binomial")
u_learners <- v_learners <- Lrnr_hal9001$new(family = "gaussian")

# simulate smaller data set for computing estimates
n_samp <- 5000
data <- sim_medoutcon_data(n_obs = n_samp)
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")
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
test_that("MSE of propensity score g estimates is sufficiently low", {
  g_mse <- mean((g_out$treat_est_train$treat_pred_A_star - g(astar, w))^2)
  expect_lt(g_mse, 0.001)
})

## fit propensity score conditioning on mediators
e_out <- fit_treat_mech(
  train_data = data, valid_data = data, contrast = c(aprime, astar),
  learners = e_learners, w_names = w_names,
  m_names = m_names, type = "e"
)
test_that("Estimates of mediator propensity score are close to the truth", {
  expect_equal(e_out$treat_est_train$treat_pred_A_prime, e(aprime, m, w),
    tol = 0.05
  )
})
test_that("MSE of mediator propensity score estimates is sufficiently low", {
  e_mse <- mean((e_out$treat_est_train$treat_pred_A_prime -
    e(aprime, m, w))^2)
  expect_lt(e_mse, 0.001)
})

## fit outcome regression
m_out <- fit_m_mech(
  train_data = data, valid_data = data, contrast = c(aprime, astar),
  learners = m_learners, m_names = m_names, w_names = w_names
)
test_that("Estimates of outcome regression m are close to the truth", {
  expect_equal(m_out$m_est_train$m_pred_A_prime, my(m, z, aprime, w),
    tol = 0.05
  )
})
test_that("MSE of outcome regression m estimates is sufficiently low", {
  m_mse <- mean((m_out$m_est_train$m_pred_A_prime - my(m, z, aprime, w))^2)
  expect_lt(m_mse, 0.005)
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
test_that("Estimates of MOC regression q are close to the truth", {
  expect_equal(q_out$moc_est_train_Z_one$moc_pred_A_prime, pz(1, aprime, w),
    tol = 0.07
  )
})
test_that("MSE of MOC regression q estimates is sufficiently low", {
  q_mse <- mean((q_out$moc_est_train_Z_one$moc_pred_A_prime -
                 pz(1, aprime, w))^2)
  expect_lt(q_mse, 0.002)
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
test_that("Estimates of MOC regression r are close to the truth", {
  expect_equal(r_out$moc_est_train_Z_one$moc_pred_A_prime,
               r(1, aprime, m, w),
    tol = 0.06
  )
})
test_that("MSE of MOC regression r estimates is sufficiently low", {
  r_mse <- mean((r_out$moc_est_train_Z_one$moc_pred_A_prime -
                 r(1, aprime, m, w))^2)
  expect_lt(r_mse, 0.001)
})

# data for fitting nuisance parameters with pseudo-outcomes
data_a_prime <- data.table::copy(data)[, A := aprime]
data_a_star <- data.table::copy(data)[, A := astar]

# the nuisance parameter u is poorly estimated in this example
if (FALSE) {
## fit u
u_out <- fit_nuisance_u(
  train_data = data,
  valid_data = data_a_prime,
  learners = u_learners,
  m_out = m_out,
  q_out = q_out,
  r_out = r_out,
  e_out = e_out,
  g_out = g_out,
  w_names = w_names
)
test_that("Estimates of pseudo-outcome regression are close to the truth", {
  expect_equal(u_out$u_pred, u(z, w), tol = 0.05)
})
test_that("MSE of pseudo-outcome regression estimates is sufficiently low", {
  u_mse <- mean((u_out$u_pred - u(z, w))^2)
  expect_lt(u_mse, 0.003)
})
}

## fit v
v_out <- fit_nuisance_v(
  train_data = data,
  valid_data = data_a_star,
  contrast = contrast,
  learners = v_learners,
  m_out = m_out,
  q_out = q_out,
  m_names = m_names,
  w_names = w_names
)
v_star <- intv(1, w) * pmaw(1, astar, w) + intv(0, w) * pmaw(0, astar, w)
test_that("Estimates of pseudo-outcome regression are close to the truth", {
  expect_equal(v_out$v_pred, v_star, tol = 0.05)
})
test_that("MSE of pseudo-outcome regression estimates is sufficiently low", {
  v_mse <- mean((v_out$v_pred - v_star)^2)
  expect_lt(v_mse, 0.002)
})
test_that("Estimates of pseudo-outcome used in v are close to the truth", {
  expect_equal(v_out$v_pseudo, intv(m, w), tol = 0.05)
})
test_that("MSE of pseudo-outcome used in v estimation is sufficiently low", {
  v_pseudo_mse <- mean((v_out$v_pseudo - intv(m, w))^2)
  expect_lt(v_pseudo_mse, 0.002)
})

