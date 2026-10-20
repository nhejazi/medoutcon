# Confidence intervals for natural/interventional (in)direct effect estimates

Compute confidence intervals for objects of class `medoutcon`, which
contain estimates produced by
[`medoutcon`](https://codex.nimahejazi.org/medoutcon/reference/medoutcon.md).

## Usage

``` r
# S3 method for class 'medoutcon'
confint(object, parm = seq_len(object$theta), level = 0.95, ...)
```

## Arguments

- object:

  An object of class `medoutcon`, as produced by invoking
  [`medoutcon`](https://codex.nimahejazi.org/medoutcon/reference/medoutcon.md),
  for which a confidence interval is to be computed.

- parm:

  A `numeric` vector indicating indices of `object$est` for which to
  return confidence intervals.

- level:

  A `numeric` indicating the level of the confidence interval to be
  computed.

- ...:

  Other arguments. Not currently used.
