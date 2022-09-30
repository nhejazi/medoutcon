test_that("two_phase_eif returns an uncentered EIF", {

  # generate fake inputs
  S <- c(1, 0, 0, 0, 1, 1)
  two_phase_weights <- c(rep(1 / 2, 3), rep(2, 3))
  eif <- c(rep(1/2, 3), rep(-1 / 2, 3))
  eif_predictions <- c(rep(1 / 3, 3), rep(-1 / 3, 3))
  plugin_est <- 1

  # independently compute the adjusted EIF
  uncentered_eif <- S * two_phase_weights * eif +
    (1 - S * two_phase_weights) * eif_predictions - plugin_est

  # assert that result is identical to two_phase_eif's
  expect_equal(
    two_phase_eif(S = S,
                  two_phase_weights = two_phase_weights,
                  eif = eif,
                  eif_predictions = eif_predictions,
                  plugin_est = plugin_est),
    uncentered_eif
  )
})
