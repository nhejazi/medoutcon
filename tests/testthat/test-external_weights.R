context("Estimator performance does not degrade with use of external weights")
source("eif_utils.R")
source("data_utils.R")

# packages
library(data.table)
library(stringr)
library(tibble)
library(hal9001)
library(sl3)

# options
set.seed(27158)
contrast <- c(0, 1)
aprime <- contrast[1]
astar <- contrast[2]
n_obs <- 5000

# 1) set up learners for nuisance parameters
hal_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", max_degree = NULL, n_folds = 5, family = "gaussian",
  lambda.min.ratio = 1 / n_obs
)
g_learners <- e_learners <- m_learners <- q_learners <- r_learners <-
  u_learners <- v_learners <- hal_lrnr

# 2) get data and column names for sl3 tasks (for convenience)
data <- sim_medoutcon_data(n_obs = n_obs)
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")
ext_wts <- runif(n_obs)
ext_wts_norm <- ext_wts / sum(ext_wts) # should we use normalized weights?

# 3) test different estimators
theta_os <- medoutcon(
  W = data[, ..w_names], A = data$A, Z = data$Z,
  M = data[, ..m_names], Y = data$Y,
  ext_weights = ext_wts,
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

# 4) compute efficient influence function based on observed data
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
psi_os <- weighted.mean(eif, ext_wts)
eif_est <- eif * ext_wts
var_eif <- var(eif_est) / n_obs

# 6) testing
test_that("Parameter estimate close to independent EIF estimates", {
  expect_equal(theta_os$theta, psi_os, tol = 0.001)
})

test_that("Variance estimate close to independent EIF variance", {
  expect_equal(theta_os$var, var_eif, tol = 0.001)
})

test_that("Mean of estimated EIF close to that of independent EIF", {
  expect_equal(abs(mean(theta_os$eif)), abs(mean(eif_est - psi_os)),
               tol = 1e-3
  )
})
