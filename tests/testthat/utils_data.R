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
  m_names <- "M"

  ## outcome
  y <- rbinom(n_obs, 1, my(m, z, a, w))

  ## construct data for output
  dat <- as.data.table(cbind(w = w, a = a, z = z, m = m, y = y))
  setnames(dat, c(w_names, "A", "Z", m_names, "Y"))
  return(dat)
}
