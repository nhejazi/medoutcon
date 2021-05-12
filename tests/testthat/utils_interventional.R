# convenience functions for nuisance parameters
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

r <- function(z, a, m, w) {
  pm(m, z, a, w) * pz(z, a, w) / pmaw(m, a, w)
}

e <- function(a, m, w) {
  pmaw(m, a, w) * g(a, w) / pmw(m, w)
}

my <- function(m, z, a, w) {
  plogis(1 / (rowSums(w) - z + a + m))
}

# compute nuisance functions by their definitions
u <- function(z, w, aprime, astar) {
  my(1, z, aprime, w) * pmaw(1, astar, w) +
    my(0, z, aprime, w) * pmaw(0, astar, w)
}

intu <- function(w, aprime, astar) {
  u(1, w, aprime, astar) * pz(1, aprime, w) + u(0, w, aprime, astar) *
    pz(0, aprime, w)
}

intv <- function(m, w, aprime) {
  my(m, 1, aprime, w) * pz(1, aprime, w) + my(m, 0, aprime, w) *
    pz(0, aprime, w)
}

sim_medoutcon_data <- function(n_obs = 1000) {
  ## baseline covariate -- simple, binary
  w_1 <- rbinom(n_obs, 1, 0.6)
  w_2 <- rbinom(n_obs, 1, 0.3)
  w_3 <- rbinom(n_obs, 1, pmin(0.2 + (w_1 + w_2) / 3, 1))
  w <- cbind(w_1, w_2, w_3)
  w_names <- paste("W", seq_len(ncol(w)), sep = "_")

  ## exposure/treatment
  a <- rbinom(n_obs, 1, g(1, w))

  ## mediator-outcome confounder affected by treatment
  z <- rbinom(n_obs, 1, pz(1, a, w))

  ## mediator (possibly multivariate)
  m <- rbinom(n_obs, 1, pm(1, z, a, w))
  m_names <- "M"

  ## outcome
  y <- rbinom(n_obs, 1, my(m, z, a, w))

  ## construct data for output
  dat <- as.data.table(cbind(w = w, a = a, z = z, m = m, y = y))
  setnames(dat, c(w_names, "A", "Z", m_names, "Y"))
  return(dat)
}
