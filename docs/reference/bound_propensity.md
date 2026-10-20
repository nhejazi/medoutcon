# Bounding propensity scores

Bounds estimated propensity score values to be within a specified range.

## Usage

``` r
bound_propensity(vals, bounds = c(0.001, 0.999))
```

## Arguments

- vals:

  A `numeric` vector of values in the interval \[0, 1\].

- bounds:

  A `numeric` vector containing two values, the first being the minimum
  allowable value and the second being the maximum allowable for values
  appearing in the vector `vals` (the previous argument).
