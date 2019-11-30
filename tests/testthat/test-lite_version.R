context("Independent estimator implementations match (package v. lite)")
source("eif_utils.R")
source("data_utils.R")

# packages
library(data.table)
library(stringr)
library(tibble)
library(hal9001)
library(SuperLearner)
library(sl3)

# options
set.seed(7128816)
contrast <- c(0, 1)
aprime <- contrast[1]
astar <- contrast[2]
n_obs <- 1000

# 1) set up learners for nuisance parameters
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
theta_os_de <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = data$Z,
  M = data[, ..m_names], Y = data$Y,
   effect = "direct",
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

theta_os_ie <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = data$Z,
  M = data[, ..m_names], Y = data$Y,
   effect = "indirect",
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

theta_tmle_de <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = data$Z,
  M = data[, ..m_names], Y = data$Y,
   effect = "direct",
  g_learners = g_learners,
  e_learners = e_learners,
  m_learners = m_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  estimator = "tmle",
  estimator_args = list(cv_folds = 2)
)

theta_tmle_ie <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = data$Z,
  M = data[, ..m_names], Y = data$Y,
   effect = "indirect",
  g_learners = g_learners,
  e_learners = e_learners,
  m_learners = m_learners,
  q_learners = q_learners,
  r_learners = r_learners,
  u_learners = u_learners,
  v_learners = v_learners,
  estimator = "tmle",
  estimator_args = list(cv_folds = 2)
)

## now fit estimators via lite version of package
candidates <- "SL.hal9001"
lite <- mediation(data,
                  weights = rep(1, nrow(data)),
                  candidatesg = candidates,
                  candidatese = candidates,
                  candidatesm = candidates,
                  candidatesr = candidates,
                  candidatesq = candidates,
                  candidatesu = candidates,
                  candidatesv = candidates,
                  nfolds = 2,
                  family.outcome = binomial())
lite


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

