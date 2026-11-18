# Marginal Coefficients from Generalized Linear Mixed Models

Calculates marginal coefficients and their standard errors from fitted
generalized linear mixed models.

## Usage

``` r
marginal_coefs(object, ...)

# S3 method for class 'MixMod'
marginal_coefs(object, std_errors = FALSE, 
   link_fun = NULL, M = 3000, K = 100, seed = 1, 
   cores = max(parallel::detectCores() - 1, 1), 
   sandwich = FALSE, ...)
```

## Arguments

- object:

  an object inheriting from class `"MixMod"`.

- std_errors:

  logical indicating whether standard errors are to be computed.

- link_fun:

  a function transforming the mean of the repeated measurements outcome
  to the linear predictor scale. Typically, this derived from the
  `family` argument of
  [`mixed_model`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.md).

- M:

  numeric scalar denoting the number of Monte Carlo samples.

- K:

  numeric scalar denoting the number of samples from the sampling
  distribution of the maximum likelihood estimates.

- seed:

  integer denoting the seed for the random number generation.

- cores:

  integer giving the number of cores to use; applicable only when
  `std_errors = TRUE`.

- sandwich:

  logical; if `TRUE` robust/sandwich standard errors are used in the
  calculations.

- ...:

  extra arguments; currently none is used.

## Details

It uses the approach of Hedeker et al. (2017) to calculate marginal
coefficients from mixed models with nonlinear link functions. The
marginal probabilities are calculated using Monte Carlo integration over
the random effects with `M` samples, by sampling from the estimated
prior distribution, i.e., a multivariate normal distribution with mean 0
and covariance matrix \\\hat{D}\\, where \\\hat{D}\\ denotes the
estimated covariance matrix of the random effects.

To calculate the standard errors, the Monte Carlo integration procedure
is repeated `K` times, where each time instead of the maximum likelihood
estimates of the fixed effects and the covariance matrix of the random
effects, a realization is used from the sampling distribution of the
maximum likelihood estimates. To speed-up this process, package
**parallel** is used.

## Value

A list of class `"m_coefs"` with components `betas` the marginal
coefficients, and when `std_errors = TRUE`, the extra components
`var_betas` the estimated covariance matrix of the marginal
coefficients, and `coef_table` a numeric matrix with the estimated
marginal coefficients, their standard errors and corresponding p-values
using the normal approximation.

## References

Hedeker, D., du Toit, S. H., Demirtas, H. and Gibbons, R. D. (2018), A
note on marginalization of regression parameters from mixed models of
binary outcomes. *Biometrics* **74**, 354–361. doi:10.1111/biom.12707

## Author

Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>

## See also

[`mixed_model`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.md)

## Examples

``` r
# \donttest{
# simulate some data
set.seed(123L)
n <- 500
K <- 15
t.max <- 25

betas <- c(-2.13, -0.25, 0.24, -0.05)
D <- matrix(0, 2, 2)
D[1:2, 1:2] <- c(0.48, -0.08, -0.08, 0.18)

times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max)))))
group <- sample(rep(0:1, each = n/2))
DF <- data.frame(year = times, group = factor(rep(group, each = K)))
X <- model.matrix(~ group * year, data = DF)
Z <- model.matrix(~ year, data = DF)

b <- cbind(rnorm(n, sd = sqrt(D[1, 1])), rnorm(n, sd = sqrt(D[2, 2])))
id <- rep(1:n, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ]))
DF$y <- rbinom(n * K, 1, plogis(eta.y))
DF$id <- factor(id)

################################################

fm1 <- mixed_model(fixed = y ~ year * group, random = ~ 1 | id, data = DF,
                   family = binomial())

fixef(fm1)                   
#>  (Intercept)         year       group1  year:group1 
#> -2.875293208  0.204144801 -0.648693468 -0.006044763 
marginal_coefs(fm1)
#> (Intercept)        year      group1 year:group1 
#>     -1.3589      0.0966     -0.3298     -0.0017 
marginal_coefs(fm1, std_errors = TRUE, cores = 1L)
#>             Estimate Std.Err  z-value p-value
#> (Intercept)  -1.3589  0.1132 -12.0052 < 1e-04
#> year          0.0966  0.0049  19.5229 < 1e-04
#> group1       -0.3298  0.1593  -2.0701 0.03844
#> year:group1  -0.0017  0.0060  -0.2795 0.77983
#> 
# }
```
