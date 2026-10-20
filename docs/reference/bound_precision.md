# Bounding to numerical precision

Bounds extreme values to a specified tolerance level, for use with
sensitive quantities that must be transformed, e.g., via
[`qlogis`](https://rdrr.io/r/stats/Logistic.html).

## Usage

``` r
bound_precision(vals, tol = 1e-06)
```

## Arguments

- vals:

  A `numeric` vector of values in the unit interval \[0, 1\].

- tol:

  A `numeric` indicating the tolerance limit to which extreme values
  should be truncated. Realizations of `val` less than `tol` are
  truncated to `tol` while those greater than (1 - `tol`) are truncated
  to (1 - `tol`).
