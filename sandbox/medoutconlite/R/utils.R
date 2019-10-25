#' Estimate treatment mechanism
#'
#' @param a ...
#' @param w ...
#'
#' @export
g <- function(a, w) {
    pscore <- (rowSums(w) / 4 + 0.1)
    return(a * pscore + (1 - a) * (1 - pscore))
}

#' Estimate intermediate confounding mechanism
#'
#' @param z ...
#' @param a ...
#' @param w ...
#'
#' @export
pz <- function(z, a, w) {
    prob1 <- plogis(rowMeans(-log(2) + w - a) + 0.2)
    return(z * prob1 + (1 - z) * (1 - prob1))
}


#' Estimate mediation mechanism
#'
#' @param m ...
#' @param z ...
#' @param a ...
#' @param w ...
#'
#' @export
pm <- function(m, z, a, w) {
    prob1 <- plogis(rowSums(log(3) * w[, -3] + a - z))
    return(m * prob1 + (1 - m) * (1 - prob1))
}

#' Estimate mediation under intervention on intermediate confounder
#'
#' @param m ...
#' @param a ...
#' @param w ...
#'
#' @export
pmaw <- function(m, a, w) {
  pm(m, 1, a, w) * pz(1, a, w) + pm(m, 0, a, w) * pz(0, a, w)
}

#' Estimate mediation under joint intervention on treatment and confounder
#'
#' @param m ...
#' @param w ...
#'
#' @export
pmw <- function(m, w) {
  pmaw(m, 1, w) * g(1, w) + pmaw(m, 0, w) * g(0, w)
}

#' Estimate re-parameterized confounding mechanism conditional on mediators
#'
#' @param z ...
#' @param a ...
#' @param m ...
#' @param w ...
#'
#' @export
r <- function(z, a, m, w) {
  pm(m, z, a, w) * pz(z, a, w) / pmaw(m, a, w)
}

#' Estimate treatment mechanism conditional on mediators
#'
#' @param a ...
#' @param m ...
#' @param w ...
#'
#' @export
e <- function(a, m, w) {
  pmaw(m, a, w) * g(a, w) / pmw(m, w)
}

#' Estimate outcome mechanism
#'
#' @param m ...
#' @param z ...
#' @param a ...
#' @param w ...
#'
#' @export
my <- function(m, z, a, w) {
  plogis(1 / as.numeric(unlist(rowSums(w) - z + a - m)))
}

#' Estimate nuisance function U
#'
#' @param z ...
#' @param w ...
#' @param aprime ...
#' @param astar ...
#'
#' @export
u <- function(z, w, aprime, astar) {
  my(1, z, aprime, w) * pmaw(1, astar, w) +
      my(0, z, aprime, w) * pmaw(0, astar, w)
}

#' Estimate integral of nuisance function U over Z
#'
#' @param w ...
#' @param aprime ...
#' @param astar ...
#'
#' @export
intu <- function(w, aprime, astar) {
  u(1, w, aprime, astar) * pz(1, aprime, w) + u(0, w, aprime, astar) *
      pz(0, aprime, w)
}

#' Estimate integral of nuisance function V over Z
#'
#' @param m ...
#' @param w ...
#' @param aprime ...
#'
#' @export
intv <- function(m, w, aprime) {
  my(m, 1, aprime, w) * pz(1, aprime, w) + my(m, 0, aprime, w) *
      pz(0, aprime, w)
}

#' Simulate example data
#'
#' n_obs ...
#'
#' @importFrom stats rbinom
#'
#' @export
simdata <- function(n_obs = 1000) {
  # helper functions for nuisance parameter estimation
  # baseline covariate -- simple, binary
  w_1 <- stats::rbinom(n_obs, 1, prob = 0.6)
  w_2 <- stats::rbinom(n_obs, 1, prob = 0.3)
  w_3 <- stats::rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
  w <- cbind(w_1, w_2, w_3)
  colnames(w) <- paste("W", seq_len(ncol(w)), sep = "")

  # exposure/treatment
  a <- as.numeric(stats::rbinom(n_obs, 1, prob = g(1, w)))

  # mediator-outcome confounder affected by treatment
  z <- stats::rbinom(n_obs, 1, pz(1, a, w))

  # mediator (possibly multivariate)
  m <- stats::rbinom(n_obs, 1, pm(1, z, a, w))
  m_names <- "M"

  # outcome
  y <- stats::rbinom(n_obs, 1, my(m, z, a, w))

  # construct data for output
  dat <- as.data.frame(cbind(W = w, A = a, Z = z, M = m, Y = y))
  return(dat)
}

#' Get truth under simulation setting
#'
#' n ...
#' contrast ...
#'
#' @importFrom stats var
#' @importFrom stringr str_subset
#' @importFrom tibble as_tibble
#'
#' @export
truth <- function(n, contrast) {
  # contrasts
  aprime <- contrast[1]
  astar <- contrast[2]

  # data
  data <- simdata(n_obs = n)
  w_names <- stringr::str_subset(colnames(data), "W")
  m_names <- stringr::str_subset(colnames(data), "M")
  w <- tibble::as_tibble(data[, w_names])
  a <- data$A
  z <- data$Z
  m <- tibble::as_tibble(data[, m_names])
  y <- data$Y

  ## compute quasi-substitution estimator of parameter
  v <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) *
      pmaw(0, astar, w)

  ## compute influence function with convenience functions
  eif <- (a == aprime) / g(aprime, w) * pmaw(m, astar, w) /
      pm(m, z, aprime, w) * (y - my(m, z, aprime, w)) + (a == aprime) /
      g(aprime, w) * (u(z, w, aprime, astar) - intu(w, aprime, astar)) +
      (a == astar) / g(astar, w) *
      (intv(m, w, aprime) - v) + v
  eif <- as.numeric(unlist(eif))

  ## computer parameter estimate and efficiency bound
  psi_os <- mean(eif)
  var_eif <- stats::var(eif)
  return(list(psi = psi_os, var = var_eif))
}
