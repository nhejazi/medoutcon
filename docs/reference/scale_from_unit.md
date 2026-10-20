# Scale values from the unit interval to their original scale

Scale values from the unit interval to their original scale

## Usage

``` r
scale_from_unit(scaled_vals, max_orig, min_orig)
```

## Arguments

- scaled_vals:

  A `numeric` vector of values scaled to lie in the closed interval \[0,
  1\] by use of
  [`scale_to_unit`](https://codex.nimahejazi.org/medoutcon/reference/scale_to_unit.md).

- max_orig:

  The maximum of the values on the original scale.

- min_orig:

  The minimum of the values on the original scale.
