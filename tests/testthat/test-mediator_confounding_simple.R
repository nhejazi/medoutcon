context("Estimators agree under mediator-outcome confounding")

library(data.table)
library(stringr)
library(future)
library(hal9001)
library(sl3)
set.seed(7128816)

################################################################################
# setup learners for the nuisance parameters
################################################################################

# instantiate some learners
mean_lrnr <- Lrnr_mean$new()
fglm_contin_lrnr <- Lrnr_glm_fast$new()
fglm_binary_lrnr <- Lrnr_glm_fast$new(family = binomial())
hal_contin_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5
)
hal_binary_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5,
  family = "binomial"
)

################################################################################
# setup data and simulate to test with estimators
################################################################################
sim_medoutcon_data <- function(n_obs = 1000) {
  # baseline covariate -- simple, binary
  w_1 <- rbinom(n_obs, 1, prob = 0.6)
  w_2 <- rbinom(n_obs, 1, prob = 0.3)
  w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
  w <- cbind(w_1, w_2, w_3)

  # exposure/treatment
  a <- as.numeric(rbinom(n_obs, 1, prob = (rowSums(w) / 4 + 0.1)))

  # mediator-outcome confounder affected by treatment
  z_0 <- rbinom(n_obs, 1, plogis(rowMeans(-log(2) + w - a) + 0.2))
  z_1 <- rbinom(n_obs, 1, plogis(rowMeans(log(3) - w[, -3] * 0.1 + a) / 3))
  z <- as.numeric(ifelse(a == 1, z_1, z_0))

  # mediator (possibly multivariate)
  m_0 <- rbinom(n_obs, 1, plogis(rowSums(-log(3) * w[, -3] - a + z)))
  m_1 <- rbinom(n_obs, 1, plogis(log(10) * w[, 1] + a - z - 0.1))
  m <- cbind(m_0, m_1)
  #m <- ifelse(z == 1, m_1, m_0)

  # outcome
  y <- rbinom(n_obs, 1, plogis(1 / (log(7) * rowSums(w) - z * a + rowSums(m))))

  # construct data for output
  dat <- as.data.table(cbind(w = w, a = a, z = z, m = m, y = y))
  return(dat)
}

# get data and column names for sl3 tasks (for convenience)
data <- sim_medoutcon_data()
w_names <- str_subset(colnames(data), "w")
z_names <- str_subset(colnames(data), "m")

################################################################################
# test different estimators
################################################################################

theta_os <- medshift_moc(
  W = data[, ..w_names], A = data$a, L = data$z,
  Z = data[, ..z_names], Y = data$y,
  contrast = c(0, 1),
  g_lrnrs = hal_binary_lrnr,
  e_lrnrs = hal_binary_lrnr,
  q_lrnrs = hal_binary_lrnr,
  r_lrnrs = hal_binary_lrnr,
  u_lrnrs = hal_contin_lrnr,
  v_lrnrs = hal_contin_lrnr,
  estimator = "onestep",
)
theta_os

#test_that("IPW and efficient one-step estimator agree", {
  #expect_equal(theta_ipw$theta, theta_aipw$theta, tol = 1e-2)
#})
