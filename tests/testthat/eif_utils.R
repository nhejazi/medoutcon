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
