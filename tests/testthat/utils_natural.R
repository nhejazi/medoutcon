###############################################################################
# from https://github.com/tlverse/tmle3mediate/blob/master/tests/testthat/helper-dgp.R
# data generating mechanism function factory
###############################################################################

# return generating functions for DGP
make_nide_dgp <- function() {
  # treatment mechanism
  g_mech <- function(w) {
    p_score <- (rowSums(w) / 4) + 0.1
    return(p_score)
  }

  # mediation mechanism
  z_mech <- function(w, a) {
    w <- data.frame(w)
    z1_prob <- 1 - plogis((a + w[, 1]) / (a + w[, 1]^3 + 0.5))
    z2_prob <- plogis((a - 1) + w[, 2] / (w[, 3] + 3))
    z3_prob <- plogis((a - 1) + 2 * w[, 1]^3 - 1 / (2 * w[, 1] + 0.5))
    return(list(
      z1_prob, z2_prob, z3_prob
    ))
  }

  # outcome mechanism for continuous
  m_mech_cont <- function(w, a, z, eps_sd = 0.5) {
    w <- data.frame(w)
    z <- data.frame(z)
    y <- z[, 1]^2 + z[, 2]^2 - z[, 3] + exp(a + z[, 3] / (1 + rowSums(w)^2)) +
      rnorm(length(a), mean = 0, sd = eps_sd)
    return(y)
  }

  # outcome mechanism for binary
  m_mech_binary <- function(w, a, z_probs, eps_sd = 0.5) {
    w <- data.frame(w)
    z_probs <- data.frame(z_probs)
    y_probs <- plogis(z_probs[, 1]^2 + z_probs[, 2]^2 - z_probs[, 3] +
      exp(a + z_probs[, 3] / (1 + rowSums(w)^2)) +
      rnorm(length(a), mean = 0, sd = eps_sd))
    return(y_probs)
  }

  # return DGP functions
  return(list(
    g_mech = g_mech,
    z_mech = z_mech,
    m_mech_cont = m_mech_cont,
    m_mech_binary = m_mech_binary
  ))
}

# just simulate baseline covariates
make_nide_baseline <- function(n_obs = 10000) {
  # baseline covariate -- simple, binary
  W_1 <- rbinom(n_obs, 1, prob = 0.5)
  W_2 <- rbinom(n_obs, 1, prob = 0.65)
  W_3 <- rbinom(n_obs, 1, prob = 0.35)
  W <- cbind(W_1, W_2, W_3)
  return(W)
}

# simulate observed data: O = (W, A, Z, Y)
make_nide_data <- function(n_obs = 10000, binary_outcome = FALSE) {
  # get data generating process functions
  dgp <- make_nide_dgp()

  # simulate baseline covariates
  W <- make_nide_baseline(n_obs)

  # get probabilities of treatment
  g_mech <- dgp$g_mech(W)

  # create treatment based on baseline W
  A <- as.numeric(rbinom(n_obs, 1, prob = g_mech))

  # get mediator probabilities
  z_mech <- dgp$z_mech(W, A)

  # simulate mediators
  Z_1 <- rbinom(n_obs, 1, prob = z_mech[[1]])
  Z_2 <- rbinom(n_obs, 1, prob = z_mech[[2]])
  Z_3 <- rbinom(n_obs, 1, prob = z_mech[[3]])
  Z <- cbind(Z_1, Z_2, Z_3)

  # create outcome as a function of A, Z, W + white noise
  # if y needs to be binary
  if (binary_outcome) {
    y_probs <- dgp$m_mech_binary(
      W, A,
      z_probs = cbind(z_mech[[1]], z_mech[[2]], z_mech[[3]]),
      eps_sd = 0.5
    )
    Y <- rbinom(n_obs, 1, prob = y_probs)
  } else {
    Y <- dgp$m_mech_cont(W, A, Z, eps_sd = 0.5)
    Y_test <- dgp$m_mech_cont(W, A, z_mech, eps_sd = 0.5)
  }

  # full data structure
  data <- as.data.table(cbind(Y, Z, A, W))
  setnames(data, c(
    "Y", paste("Z", 1:3, sep = "_"), "A",
    paste("W", seq_len(dim(W)[2]), sep = "_")
  ))
  return(data)
}

###############################################################################
# https://github.com/tlverse/tmle3mediate/blob/master/tests/testthat/helper-NDE_NIE.R
###############################################################################
get_truth_nide <- function(n_obs = 1e6, binary_outcome = FALSE, EIC = FALSE) {
  # get DGP generating functions
  dgp <- make_nide_dgp()

  # compute baseline data
  W <- make_nide_baseline(n_obs)

  # get treatment
  g1 <- dgp$g_mech(W)
  A <- rbinom(n_obs, 1, prob = g1)
  g <- A * g1 + (1 - A) * (1 - g1)

  # simulate mediators
  z_probs <- dgp$z_mech(W, A)
  Z_1 <- rbinom(n_obs, 1, prob = z_probs[[1]])
  Z_2 <- rbinom(n_obs, 1, prob = z_probs[[2]])
  Z_3 <- rbinom(n_obs, 1, prob = z_probs[[3]])
  Z <- cbind(Z_1, Z_2, Z_3)

  # get mediator counterfactuals
  Z_1_probs <- dgp$z_mech(W, 1)
  Z_0_probs <- dgp$z_mech(W, 0)
  Z1_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[1]])
  Z1_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[1]])
  Z2_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[2]])
  Z2_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[2]])
  Z3_1 <- rbinom(n_obs, 1, prob = Z_1_probs[[3]])
  Z3_0 <- rbinom(n_obs, 1, prob = Z_0_probs[[3]])
  Z1 <- cbind(Z_1 = Z1_1, Z_2 = Z2_1, Z_3 = Z3_1)
  Z0 <- cbind(Z_1 = Z1_0, Z_2 = Z2_0, Z_3 = Z3_0)

  # compute outcome counterfactuals: Y(1, 1), Y(0, 0), Y(1, 0), Y(0, 1)
  if (binary_outcome) {
    EY_A1_Z1 <- dgp$m_mech_binary(W, 1, Z_1_probs, eps_sd = 0)
    EY_A0_Z0 <- dgp$m_mech_binary(W, 0, Z_0_probs, eps_sd = 0)
    EY_A1_Z0 <- dgp$m_mech_binary(W, 1, Z_0_probs, eps_sd = 0)
    EY_A0_Z1 <- dgp$m_mech_binary(W, 0, Z_1_probs, eps_sd = 0)

    # compute TRUE M under counterfactual regimes
    m_Ais1 <- dgp$m_mech_binary(W, 1, z_probs, eps_sd = 0)
    m_Ais0 <- dgp$m_mech_binary(W, 0, z_probs, eps_sd = 0)

    # compute difference of counterfactuals based on substitution formula
    psi_Z_NDE <- EY_A1_Z0 - EY_A0_Z0
    psi_Z_NIE <- EY_A1_Z1 - EY_A1_Z0
  } else {
    EY_A1_Z1 <- dgp$m_mech_cont(W, 1, Z1, eps_sd = 0)
    EY_A0_Z0 <- dgp$m_mech_cont(W, 0, Z0, eps_sd = 0)
    EY_A1_Z0 <- dgp$m_mech_cont(W, 1, Z0, eps_sd = 0)
    EY_A0_Z1 <- dgp$m_mech_cont(W, 0, Z1, eps_sd = 0)

    # compute TRUE M under counterfactual regimes
    m_Ais1 <- dgp$m_mech_cont(W, 1, Z, eps_sd = 0)
    m_Ais0 <- dgp$m_mech_cont(W, 0, Z, eps_sd = 0)

    # compute difference of counterfactuals based on substitution formula
    psi_Z_NDE <- EY_A1_Z0 - EY_A0_Z0
    psi_Z_NIE <- EY_A1_Z1 - EY_A1_Z0
  }

  # compute approximate natural (in)direct effects
  NDE <- mean(psi_Z_NDE)
  NIE <- mean(psi_Z_NIE)

  if (EIC) {
    df1 <- as.data.frame(cbind(Z1, W))
    df0 <- as.data.frame(cbind(Z0, W))
    Q_Z1 <- df1 %>%
      dplyr::group_by(
        Z_1, Z_2, Z_3,
        W_1, W_2, W_3
      ) %>%
      summarize(
        prob = n() / n_obs
      )
    Q_Z0 <- df0 %>%
      dplyr::group_by(
        Z_1, Z_2, Z_3,
        W_1, W_2, W_3
      ) %>%
      summarize(
        prob = n() / n_obs
      )
    Q_Z_Ais1 <- df1 %>%
      dplyr::left_join(
        Q_Z1
      ) %>%
      pull(prob)
    Q_Z_Ais0 <- df0 %>%
      dplyr::left_join(
        Q_Z0
      ) %>%
      pull(prob)
    Q_Z_ratio <- Q_Z_Ais0 / Q_Z_Ais1

    if (binary_outcome) {
      Y <- dgp$m_mech_binary(W, A, z_probs)
      Qbar_Y <- dgp$m_mech_binary(W, A, Z, eps_sd = 0)
    } else {
      Y <- dgp$m_mech_cont(W, A, Z)
      Qbar_Y <- dgp$m_mech_cont(W, A, Z, eps_sd = 0)
    }

    # natural direct effect
    HY <- (A / g1) * (Q_Z_Ais0 / Q_Z_Ais1) - ((1 - A) / (1 - g1))
    HZ <- 1 / (1 - g1)

    # compute individual scores for DY, DA, DW
    D_Y_NDE <- HY * (Y - Qbar_Y)
    D_Z_NDE <- (1 - A) * HZ * (m_Ais1 - m_Ais0 - psi_Z_NDE)
    D_W_NDE <- psi_Z_NDE - NDE
    EIC_NDE <- D_Y_NDE + D_Z_NDE + D_W_NDE

    # natural indirect effect
    EIC_NIE <- ((A / g1) * (Y - EY_A1_Z1 - Q_Z_ratio) * (Y - m_Ais1)) -
      ((1 - A) / (1 - g1)) * (m_Ais1 - EY_A1_Z0) +
      psi_Z_NIE - NIE

    # output: true values of nuisance parameters
    return(list(
      EY_A1_Z1 = EY_A1_Z1,
      EY_A1_Z0 = EY_A1_Z0,
      EY_A0_Z1 = EY_A0_Z1,
      EY_A0_Z0 = EY_A0_Z0,
      NDE = NDE,
      NIE = NIE,
      EIC_NDE = EIC_NDE,
      EIC_NIE = EIC_NIE
    ))
  } else {
    # output: true values of nuisance parameters
    return(list(
      EY_A1_Z1 = EY_A1_Z1,
      EY_A1_Z0 = EY_A1_Z0,
      EY_A0_Z1 = EY_A0_Z1,
      EY_A0_Z0 = EY_A0_Z0,
      NDE = NDE,
      NIE = NIE
    ))
  }
}
