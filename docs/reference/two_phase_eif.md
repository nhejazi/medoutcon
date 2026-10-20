# Two-phase sampling adjusted, un-centered efficient influence function

Adjust the efficient influence function to account for the use of
two-phase sampling designs in measuring the mediators.

## Usage

``` r
two_phase_eif(R, two_phase_weights, eif, eif_predictions, plugin_est)
```

## Arguments

- R:

  A `logical` vector indicating whether an sampled observation's
  mediators were measured using a two-phase sampling design.

- two_phase_weights:

  A `numeric` vector of known observation-level weights corresponding to
  the inverse probability of the mediators being measured. These weights
  should only be provided if the two-phase sampling indicator `R` is
  specified.

- eif:

  A `numeric` vector of the efficient influence function.

- eif_predictions:

  A `numeric` vector of the predicted efficient influence function,
  conditioning on the mediator being measured.

- plugin_est:

  A `numeric` corresponding to the plug-in estimate for the given
  contrast.

## Value

An un-centered efficeint influence function that accounts for the
two-phase sampling design.
