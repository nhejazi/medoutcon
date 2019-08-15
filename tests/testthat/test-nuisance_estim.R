library(tidyverse)
library(data.table)
library(hal9001)
library(sl3)
g_learners <- Lrnr_glm_fast$new(family = stats::binomial())
e_learners <- m_learners <- q_learners <- r_learners <-
  Lrnr_hal9001$new(family = "binomial")
u_learners <- v_learners <- Lrnr_hal9001$new(family = "gaussian")

sim_medoutcon_data <- function(n_obs = 1000) {
  ## baseline covariate -- simple, binary
  w_1 <- rbinom(n_obs, 1, prob = 0.6)
  w_2 <- rbinom(n_obs, 1, prob = 0.3)
  w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
  w <- cbind(w_1, w_2, w_3)
  w_names <- paste("W", seq_len(ncol(w)), sep = "_")

  ## exposure/treatment
  a <- as.numeric(rbinom(n_obs, 1, prob = g(1, w)))

  ## mediator-outcome confounder affected by treatment
  z <- rbinom(n_obs, 1, pz(1, a, w))

  ## mediator (possibly multivariate)
  m <- rbinom(n_obs, 1, pm(1, z, a, w))
  ## m_2_prob <- plogis(log(10) * w[, 1] + a - z - 0.1)
  ## m_2 <- runif(n_obs, min(m_2_prob), max(m_2_prob))
  ## m <- cbind(m_1, m_2)
  ## m_names <- paste("M", seq_len(ncol(m)), sep = "_")
  m_names <- "M"

  ## outcome
  y <- rbinom(n_obs, 1, my(m, z, a, w))

  ## construct data for output
  dat <- as.data.table(cbind(w = w, a = a, z = z, m = m, y = y))
  setnames(dat, c(w_names, "A", "Z", m_names, "Y"))
  return(dat)
}

################################################################################
# convenience functions for nuisance parameters
################################################################################
g <- function(a, w) {
  pscore <- (rowSums(w) / 4 + 0.1)
  return(a * pscore + (1 - a) * (1 - pscore))
}

pz <- function(z, a, w) {
  prob1 <- plogis(rowMeans(-log(2) + w - a) + 0.2)
  return(z * prob1 + (1 - z) * (1 - prob1))
}

pm <- function(m, z, a, w) {
  prob1 <- plogis(rowSums(log(3) * w[, -3] + a - z))
  return(m * prob1 + (1 - m) * (1 - prob1))
}

pmaw <- function(m, a, w) {
  pm(m, 1, a, w) * pz(1, a, w) +
    pm(m, 0, a, w) * pz(0, a, w)
}

pmw <- function(m, w) {
  pmaw(m, 1, w) * g(1, w) +
    pmaw(m, 0, w) * g(0, w)
}

r <- function(z, a, m, w) pm(m, z, a, w) * pz(z, a, w) / pmaw(m, a, w)

e <- function(a, m, w) pmaw(m, a, w) * g(a, w) / pmw(m, w)

my <- function(m, z, a, w) plogis(1 / (rowSums(w) - z + a - m))

# compute nuisance functions by their definitions
u <- function(z, w) {
  my(1, z, aprime, w) * pmaw(1, astar, w) +
    my(0, z, aprime, w) * pmaw(0, astar, w)
}

intu <- function(w) {
  u(1, w) * pz(1, aprime, w) + u(0, w) * pz(0, aprime, w)
}

intv <- function(m, w) {
  my(m, 1, aprime, w) * pz(1, aprime, w) + my(m, 0, aprime, w) *
    pz(0, aprime, w)
}

################################################################################
# simulate data
################################################################################

# data and helper variables
n_samp <- 1e7
data <- sim_medoutcon_data(n_obs = n_samp)
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")
w <- as_tibble(data)[, w_names]
a <- data$A
z <- data$Z
m <- data$M
y <- data$Y

# contrasts
contrast <- c(1, 0)
aprime <- contrast[1]
astar <- contrast[2]

# compute parameter estimate and influence function with convenience functions
v <- intv(1, w) * pmaw(1, astar, w) + intv(0, w) * pmaw(0, astar, w)

eif <- (a == aprime) / g(aprime, w) * pmaw(m, astar, w) /
  pm(m, z, aprime, w) * (y - my(m, z, aprime, w)) + (a == aprime) /
    g(aprime, w) * (u(z, w) - intu(w)) + (a == astar) / g(astar, w) *
    (intv(m, w) - v) + v
psi_os <- mean(eif)
var_eif <- var(eif)

################################################################################
# nuisance components for one-step estimator
################################################################################

n_samp <- 1e4
data <- sim_medoutcon_data(n_obs = n_samp)
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")
w <- as_tibble(data)[, w_names]
a <- data$A
z <- data$Z
m <- data$M
y <- data$Y

if (FALSE) {

  ## fit propensity score
  g_out <- fit_treat_mech(
    train_data = data, valid_data = data, contrast = c(aprime, astar),
    learners = g_learners, w_names = w_names, type = "g"
  )
  test_that("Estimates of propensity score are close to the truth", {
    expect_equal(g_out$treat_est_train$treat_pred_A_star, g(astar, w),
                 tol = 5e-2)
  })
  test_that("MSE of propensity score g estimates is sufficiently low", {
    g_mse <- mean((g_out$treat_est_train$treat_pred_A_star - g(astar, w))^2)
    expect_lt(g_mse, 10/n_samp)
  })

  ## fit propensity score conditioning on mediators
  e_out <- fit_treat_mech(
    train_data = data, valid_data = data, contrast = c(aprime, astar),
    learners = e_learners, w_names = w_names,
    m_names = m_names, type = "e"
  )
  test_that("Estimates of mediator propensity score are close to the truth", {
    expect_equal(e_out$treat_est_train$treat_pred_A_prime, e(aprime, m, w),
                 tol = 5e-2)
  })
  test_that("MSE of mediator propensity score estimates is sufficiently low", {
    e_mse <- mean((e_out$treat_est_train$treat_pred_A_prime -
                   e(aprime, m, w))^2)
    expect_lt(e_mse, 10 / n_samp)
  })

  ## fit outcome regression
  m_out <- fit_m_mech(
    train_data = data, valid_data = data, contrast = c(aprime, astar),
    learners = m_learners, m_names = m_names, w_names = w_names
  )
  test_that("Estimates of outcome regression m are close to the truth", {
    expect_equal(m_out$m_est_train$m_pred_A_prime, my(m, z, aprime, w),
      tol = 5e-2
    )
  })
  test_that("MSE of outcome regression m estimates is sufficiently low", {
    m_mse <- mean((m_out$m_est_train$m_pred_A_prime - my(m, z, aprime, w))^2)
    expect_lt(m_mse, 1 / n_samp)
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
    expect_equal(q_out$moc_est_train$moc_pred_A_prime, pz(1, aprime, w),
      tol = 5e-2
    )
  })
  test_that("MSE of MOC regression q estimates is sufficiently low", {
    q_mse <- mean((q_out$moc_est_train$moc_pred_A_prime - pz(1, aprime, w))^2)
    expect_lt(q_mse, 1 / n_samp)
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
    expect_equal(r_out$moc_est_train$moc_pred_A_prime, r(1, aprime, m, w),
      tol = 5e-2
    )
  })
  test_that("MSE of MOC regression r estimates is sufficiently low", {
    r_mse <- mean((r_out$moc_est_train$moc_pred_A_prime - r(1, aprime, m, w))^2)
    expect_lt(r_mse, 10 / n_samp)
  })

  ## fit u
  data_a_prime <- data.table::copy(data)[, A := aprime]
  data_a_star <- data.table::copy(data)[, A := astar]
  u_out <- fit_nuisance_u(
    train_data = data,
    valid_data = data_a_prime,
    learners = u_learners,
    m_out = m_out,
    q_out = q_out,
    r_out = r_out,
    e_out = e_out,
    w_names = w_names
  )
  u_prime <- u_out$u_pred
  test_that("Estimates of pseudo-outcome regression are close to the truth", {
    expect_equal(u_prime, u(z, w), tol = 5e-2)
  })
  test_that("MSE of pseudo-outcome regression estimates is sufficiently low", {
    u_mse <- mean((u_prime - u(z, w))^2)
    expect_lt(u_mse, 10 / n_samp)
  })

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
  test_that("Estimates of pseudo-outcome regression are close to the truth", {
    v_star <- intv(1, w) * pmaw(1, astar, w) + intv(0, w) * pmaw(0, astar, w)
    expect_equal(v_out$v_pred, v_star, tol = 5e-2)
  })
  test_that("MSE of pseudo-outcome regression estimates is sufficiently low", {
    v_mse <- mean((v_out$v_pred - v_star)^2)
    expect_lt(v_mse, 10 / n_samp)
  })
  test_that("Estimates of pseudo-outcome used in v are close to the truth", {
    expect_equal(v_out$v_pseudo, intv(m, w), tol = 5e-2)
  })
  test_that("MSE of pseudo-outcome used in v estimation is sufficiently low", {
    v_pseudo_mse <- mean((v_out$v_pseudo - intv(m, w))^2)
    expect_lt(v_pseudo_mse, 1 / n_samp)
  })

  ## extract pre-computed components from nuisance function estimates
  m_prime <- m_out$m_est_valid$m_pred_A_prime
  q_prime <- q_out$moc_est_valid$moc_pred_A_prime
  r_prime <- r_out$moc_est_valid$moc_pred_A_prime
  e_prime <- e_out$treat_est_valid$treat_pred_A_prime
  e_star <- e_out$treat_est_valid$treat_pred_A_star
  g_star <- g_out$treat_est_valid$treat_pred_A_star

  ## ad-hoc function to compute integral over Z for nuisance parameter u
  u_int_eif <- lapply(unique(data$Z), function(confounder_val) {
    # intervene on training and validation data sets
    data_z_interv <- data.table::copy(data)
    data_z_interv[, `:=`(
      Z = confounder_val,
      A = contrast[1],
      U_pseudo = u_prime
    )]

    # predict u(z, a', w) using intervened data with treatment set A = a'
    u_task_z_interv <- sl3::sl3_Task$new(
      data = data_z_interv,
      covariates = c(w_names, "A", "Z"),
      outcome_type = "continuous",
      outcome = "U_pseudo"
    )
    u_prime_z_interv <- u_out[["u_fit"]]$predict(u_task_z_interv)

    # task for predicting from trained q regression model
    q_reg_v_subtask <- sl3::sl3_Task$new(
      data = data,
      covariates = c("A", w_names),
      outcome_type = "binomial",
      outcome = "Z"
    )

    # q nuisance regression after intervening on mediator-outcome confounder
    # NOTE: for binary Z, this returns P(Z = 1 | ...) by definition but what we
    #       want is actually P(Z = z | ...)
    q_pred_z_natural <- (confounder_val * q_prime) +
      ((1 - confounder_val) * (1 - q_prime))

    # return partial pseudo-outcome for v nuisance regression
    out <- u_prime_z_interv * q_pred_z_natural
    return(out)
  })
  u_int_eif <- do.call(`+`, u_int_eif)
  test_that("Estimates of u-q integral (over z) are close to the truth", {
    expect_equal(u_int_eif, intu(w), tol = 5e-2)
  })
  test_that("MSE of u-q integral (over z) is sufficiently low", {
    intu_mse <- mean((u_int_eif - intu(w))^2)
    expect_lt(intu_mse, 1 / n_samp)
  })

  ## finally, compute un-centered influence function for one-step
  ipw_a_prime <- as.numeric(data$A == aprime) / g_star
  ipw_a_star <- as.numeric(data$A == astar) / g_star
  eif_resid_y <- (q_prime / r_prime) * (e_star / e_prime) * (data$Y - m_prime)
  eif_test <- (ipw_a_prime * eif_resid_y) +
    (ipw_a_prime * (u_prime - u_int_eif)) +
    (ipw_a_star * (v_out$v_pseudo - v_out$v_pred)) + v_out$v_pred
  psi_os_test <- mean(eif_test)
  var_eif_test <- var(eif_test)
  test_that("Estimates based on influence function are close to the truth", {
    expect_equal(eif_test, eif, tol = 5e-2)
    expect_equal(psi_os_test, psi_os, tol = 5e-2)
    expect_equal(var_eif_test, var_eif, tol = 5e-2)
  })
  test_that("MSE of estimated influence function  is sufficiently low", {
    eif_mse <- mean((eif_test - eif)^2)
    expect_lt(eif_mse, 1 / n_samp)
  })
}
