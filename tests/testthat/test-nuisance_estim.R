library(tidyverse)
library(data.table)
library(sl3)
g_learners <- q_learners <- r_learners <-
  Lrnr_glm_fast$new(family = stats::binomial())
e_learners <- m_learners <- Lrnr_hal9001$new(family = "binomial")
u_learners <- v_learners <- Lrnr_glm_fast$new()

g <- function(a, w) {
    pscore <- (rowSums(w) / 4 + 0.1)
    return(a * pscore + (1 - a) * (1 - pscore))
}

pz <- function(z, a, w) {
    prob1 <- plogis(rowMeans(-log(2) + w - a) + 0.2)
    return(z * prob1 + (1 - z) * (1 - prob1))
}

pm <- function(m, z, a, w) {
    prob1 <- plogis(rowSums(log(3) * w[,-3] + a - z))
    return(m * prob1 + (1 - m) * (1 - prob1))
}

pmaw <- function(m, a, w) {
    pm(m, 1, a, w) * pz(1, a, w) +
        pm(m, 0, a, w) * pz(0, a, w)
}

pmw <-  function(m, w) {
    pmaw(m, 1, w) * g(1, w) +
        pmaw(m, 0, w) * g(0, w)
}

r <- function(z, a, m, w) pm(m, z, a, w) * pz(z, a, w) / pmaw(m, a, w)

e <- function(a, m, w) pmaw(m, a, w) * g(a, w) / pmw(m, w)

my <- function(m, z, a, w) plogis(1 / (rowSums(w) - z + a - m))

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
    ##m_2_prob <- plogis(log(10) * w[, 1] + a - z - 0.1)
    ##m_2 <- runif(n_obs, min(m_2_prob), max(m_2_prob))
    ##m <- cbind(m_1, m_2)
    ##m_names <- paste("M", seq_len(ncol(m)), sep = "_")
    m_names <- "M"

    ## outcome
    y <- rbinom(n_obs, 1, my(m, z, a, w))

    ## construct data for output
    dat <- as.data.table(cbind(w = w, a = a, z = z, m = m, y = y))
    setnames(dat, c(w_names, "A", "Z", m_names, "Y"))
    return(dat)
}


# data and helper variables
data <- sim_medoutcon_data(n_obs = 1e3)
w_names <- str_subset(colnames(data), "W")
m_names <- str_subset(colnames(data), "M")
w <- data[, ..w_names]
a <- data$A
z <- data$Z
m <- data$M
y <- data$Y

# constrasts
contrast <- c(0, 1)
aprime <- contrast[1]
astar  <- contrast[2]

# compute nuisance functions by their definitions
u <- function(z, w) {
    my(1, z, aprime, w) * pmaw(1, astar, w) +
        my(0, z, aprime, w) * pmaw(0, astar, w)
}

intu <- function(w) {
    u(1, w) * pz(1, aprime, w) + u(0, w) * pz(0, aprime, w)
}

intv <- function(m, w){
    my(m, 1, aprime, w) * pz(1, aprime, w) +  my(m, 0, aprime, w) *
        pz(0, aprime, w)
}

v <- intv(1, w) * pmaw(1, astar, w) + intv(0, w) * pmaw(0, astar, w)

eif <- (a == aprime) / g(astar, w) * pz(z, aprime, w) / r(z, aprime, m, w) *
    e(astar, m, w) / e(aprime, m, w) * (y - my(m, z, aprime, w)) +
    (a == aprime) / g(astar, w) * (u(z, w) - intu(w)) +
    (a == astar) / g(astar, w) * (intv(m, w) - v) + v



# nuisance components for one-step estimator
g_out <- fit_g_mech(train_data = data, contrast = c(0, 1),
                    learners = g_learners, w_names = w_names)
test_that("", {
  expect_equal(g_out$g_est$g_pred_A_star, g(astar, w), tol = 5e-2)
})

e_out <- fit_e_mech(train_data = data, contrast = c(0, 1),
                    learners = e_learners, w_names = w_names,
                    m_names = m_names)
test_that("", {
  expect_equal(e_out$e_est$e_pred_A_prime, e(aprime, m, w), tol = 5e-2)
})

m_out <- fit_m_mech(train_data = data, contrast = c(0, 1),
                    learners = m_learners, m_names = m_names, w_names = w_names)
test_that("", {
  expect_equal(m_out$m_est$m_pred_A_prime, my(m, z, aprime, w), tol = 5e-2)
})
