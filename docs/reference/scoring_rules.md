# Proper Scoring Rules for Categorical Data

Calculates the logarithmic, quadratic/Brier and spherical based on a
fitted mixed model for categorical data.

## Usage

``` r
scoring_rules(object, newdata, newdata2 = NULL, max_count = 2000, 
    return_newdata = FALSE)
```

## Arguments

- object:

  an object inheriting from class `"MixMod"`.

- newdata:

  a data.frame based on which to estimate the random effect and
  calculate predictions. It should contain the response variable.

- newdata2:

  a data.frame based on which to estimate the random effect and
  calculate predictions. It should contain the response variable.

- max_count:

  numeric scalar denoting the maximum count up to which to calculate
  probabilities; this is relevant for count response data.

- return_newdata:

  logical; if `TRUE` the values of the scoring rules are ruturned as
  extra columns of the `newdata` or `newdata2` data.frame.

## Value

A data.frame with (extra) columns the values of the logarithmic,
quadratic and spherical scoring rules calculated based on the fitted
model and the observed responses in `newdata` or `newdata2`.

## References

Carvalho, A. (2016). An overview of applications of proper scoring
rules. *Decision Analysis* **13**, 223–242. doi:10.1287/deca.2016.0337

## Author

Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>

## See also

[`mixed_model`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.md),
[`predict.MixMod`](https://drizopoulos.github.io/GLMMadaptive/reference/methods.md)

## Examples

``` r
# \donttest{
NA
#> [1] NA
NA
#> [1] NA
NA
#> [1] NA
# }
```
