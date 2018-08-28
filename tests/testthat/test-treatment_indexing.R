context("Indexing by treatment vector matches multiplicative indexing")

library(data.table)
library(stringr)
set.seed(429153)
delta <- 0.5

################################################################################
# setup data and simulate to test with estimators
################################################################################

# simulate simple data for simple simulation
make_simulated_data <- function(n_obs = 1000, # number of observations
                                n_w = 3, # number of baseline covariates
                                p_w = 0.5, # prob. of success in baseline vars
                                delta_shift = delta # posited shift parameter
) {
  # baseline covariate -- simple, binary
  W <- as.matrix(replicate(n_w, rbinom(n_obs, 1, prob = p_w)))

  # create treatment based on baseline W
  A <- as.numeric(rbinom(n_obs, 1, prob = (rowSums(W) / 4 + 0.1)))

  # mediators to affect the outcome
  ## 1st mediator (binary)
  z1_prob <- (A + W[, 1]^2 + runif(n_obs, 0, 0.1)) / max(A + W[, 3] + 0.1)
  Z_1 <- rbinom(n_obs, 1, prob = z1_prob)
  ## 2nd mediator (binary)
  z2_form <- (A - 1)^3 + W[, 2] - W[, 3] + runif(n_obs, 0, 0.2)
  z2_prob <- abs(z2_form) / max(abs(z2_form))
  z2_prob <- z2_prob + runif(n_obs, 0, 1 - max(z2_prob))
  Z_2 <- rbinom(n_obs, 1, prob = z2_prob)
  ## 3rd mediator (binary)
  z3_form <- (A - 1)^2 + 2 * W[, 1]^3 + rnorm(n_obs)
  z3_prob <- abs(z3_form) / max(abs(z3_form))
  Z_3 <- rbinom(n_obs, 1, prob = z3_prob)
  ## build matrix of mediators
  Z <- cbind(Z_1, Z_2, Z_3)

  # create outcome as a linear function of A, W + white noise
  Y <- Z_1 + Z_2 - Z_3 + A - 0.1 * rowSums(W)^2 + rnorm(n_obs, mean = 0, sd = 1)

  # full data structure
  data <- as.data.table(cbind(Y, Z, A, W))
  setnames(data, c(
    "Y", paste("Z", 1:3, sep = "_"), "A",
    paste("W", seq_len(dim(W)[2]), sep = "_")
  ))
  return(data)
}

data <- make_simulated_data()
head(data)

# column names for sl3 tasks (for convenience)
z_names <- colnames(data)[str_detect(colnames(data), "Z")]
w_names <- colnames(data)[str_detect(colnames(data), "W")]


# get indices of treatment and control units
idx_treat <- which(data$A == 1)
idx_cntrl <- which(data$A == 0)


# just try indexing by the treatment mechanism
g_out <- fit_g_mech(
  data = data,
  delta = delta,
  lrnr_stack = Lrnr_hal9001$new(fit_type = "glmnet", family = "binomial"),
  w_names = w_names
)

g_shifted_treat <- g_out$g_est$g_pred_A1
g_shifted_cntrl <- g_out$g_est$g_pred_A0

# construct shifted vector based on indexing
g_shifted_idx <- rep(NA, length(g_shifted_treat))
g_shifted_idx[idx_treat] <- g_shifted_treat[idx_treat]
g_shifted_idx[idx_cntrl] <- g_shifted_cntrl[idx_cntrl]

# construct shifted vector based on treatment vector  multiplication
g_shifted_mult <- (g_shifted_treat * data$A) + g_shifted_cntrl * (1 - data$A)

# test that indexing approaches produce identical results
test_that("Different indexing approaches produce identical results", {
  expect_identical(g_shifted_idx, g_shifted_mult)
})
