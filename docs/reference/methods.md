# Various Methods for Standard Generics

Methods for object of class `"MixMod"` for standard generic functions.

## Usage

``` r
coef(object, ...)

# S3 method for class 'MixMod'
coef(object, sub_model = c("main", "zero_part"), 
    ...)

fixef(object, ...)

# S3 method for class 'MixMod'
fixef(object, sub_model = c("main", "zero_part"), ...)

ranef(object, ...)

# S3 method for class 'MixMod'
ranef(object, post_vars = FALSE, ...)

confint(object, parm, level = 0.95, ...)

# S3 method for class 'MixMod'
confint(object, 
  parm = c("fixed-effects", "var-cov","extra", "zero_part"), 
  level = 0.95, sandwich = FALSE, ...)

anova(object, ...)

# S3 method for class 'MixMod'
anova(object, object2, test = TRUE, 
  L = NULL, sandwich = FALSE, ...)

fitted(object, ...)

# S3 method for class 'MixMod'
fitted(object, 
  type = c("mean_subject", "subject_specific", "marginal"),
  link_fun = NULL, ...)

residuals(object, ...)

# S3 method for class 'MixMod'
residuals(object, 
  type = c("mean_subject", "subject_specific", "marginal"), 
  link_fun = NULL, tasnf_y = function (x) x, ...)
  
predict(object, ...)

# S3 method for class 'MixMod'
predict(object, newdata, newdata2 = NULL, 
    type_pred = c("response", "link"),
    type = c("mean_subject", "subject_specific", "marginal", "zero_part"),
    se.fit = FALSE, M = 300, df = 10, scale = 0.3, level = 0.95, 
    seed = 1, return_newdata = FALSE, sandwich = FALSE, ...)
    
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'MixMod'
simulate(object, nsim = 1, seed = NULL, 
    type = c("subject_specific", "mean_subject"), new_RE = FALSE,
    acount_MLEs_var = FALSE, sim_fun = NULL, 
    sandwich = FALSE, ...)
    
terms(x, ...)

# S3 method for class 'MixMod'
terms(x, type = c("fixed", "random", "zi_fixed", "zi_random"), ...)

formula(x, ...)

# S3 method for class 'MixMod'
formula(x, type = c("fixed", "random", "zi_fixed", "zi_random"), ...)


model.frame(formula, ...)

# S3 method for class 'MixMod'
model.frame(formula, type = c("fixed", "random", "zi_fixed", "zi_random"), ...)

model.matrix(object, ...)

# S3 method for class 'MixMod'
model.matrix(object, type = c("fixed", "random", "zi_fixed", "zi_random"), ...)

nobs(object, ...)

# S3 method for class 'MixMod'
nobs(object, level = 1, ...)

VIF(object, ...)

# S3 method for class 'MixMod'
VIF(object, type = c("fixed", "zi_fixed"), ...)

cooks.distance(model, ...)

# S3 method for class 'MixMod'
cooks.distance(model, cores = max(parallel::detectCores() - 1, 1), ...)
```

## Arguments

- object, object2, x, formula, model:

  objects inheriting from class `"MixMod"`. When `object2` is also
  provided, then the model behind `object` must be nested within the
  model behind `object2`.

- sub_model:

  character string indicating for which sub-model to extract the
  estimated coefficients; it is only relevant for zero-inflated models.

- post_vars:

  logical; if `TRUE` the posterior variances of the random effects are
  returned as an extra attribute of the numeric matrix produced by
  `ranef()`.

- parm:

  character string; for which type of parameters to calculate confidence
  intervals. Option `"var-cov"` corresponds to the variance-covariance
  matrix of the random effects. Option `extra` corresponds to extra
  (shape/dispersion) parameters in the distribution of the outcome
  (e.g., the \\\theta\\ parameter in the negative binomial family).
  Option `zero_inflated` corresponds to the coefficients of the
  zero-inflated sub-model.

- level:

  numeric scalar between 0 and 1 denoting the level of the confidence
  interval. In the `nobs()` method it denotes the level at which the
  number of observations is counted. The value 0 corresponds to the
  number of independent sample units determined by the number of levels
  of the grouping variable. If set to a value greater than zero, it
  returns the total number of observations.

- test:

  logical; should a p-value be calculated.

- L:

  a numeric matrix representing a contrasts matrix. This is only used
  when in `anova()` only `object` is provided, and it can only be
  specified for the fixed effects. When `L` is used, a Wald test is
  performed.

- sandwich:

  logical; if `TRUE` the sandwich estimator is used in the calculation
  of standard errors.

- type:

  character string indicating the type of fitted values / residuals /
  predictions / variance inflation factors to calculate. Option
  `"mean_subject"` corresponds to only using the fixed-effects part;
  option `"subject_specific"` corresponds to using both the fixed- and
  random-effects parts; option `"marginal"` is based in multiplying the
  fixed effects design matrix with the marginal coefficients obtained by
  [`marginal_coefs`](https://drizopoulos.github.io/GLMMadaptive/reference/marginal_coefs.md).

- link_fun:

  the `link_fun` of
  [`marginal_coefs`](https://drizopoulos.github.io/GLMMadaptive/reference/marginal_coefs.md).

- tasnf_y:

  a function to transform the grouped / repeated measurements outcome
  before calculating the residuals; for example, relevant in two-part
  models for semi-continuous data, in which it is assumed that the log
  outcome follows a normal distribution.

- newdata, newdata2:

  a data frame based on which predictions are to be calculated.
  `newdata2` is only relevant when `level = "subject_specific"`; see
  **Details** for more information.

- type_pred:

  character string indicating at which scale to calculate predictions.
  Options are `"link"` indicating to calculate predictions at the link
  function / linear predictor scale, and `"response"` indicating to
  calculate predictions at the scale of the response variable.

- se.fit:

  logical, if `TRUE` standard errors of predictions are returned.

- M:

  numeric scalar denoting the number of Monte Carlo samples; see
  **Details** for more information.

- df:

  numeric scalar denoting the degrees of freedom for the Student's t
  proposal distribution; see **Details** for more information.

- scale:

  numeric scalar or vector denoting the scaling applied to the
  subject-specific covariance matrices of the random effects; see
  **Details** for more information.

- seed:

  numerical scalar giving the seed to be used in the Monte Carlo scheme.

- return_newdata:

  logical; if `TRUE` the `predict()` method returns a copy of the
  `newdata` and of `newdata2` if the corresponding argument was not
  `NULL`, with extra columns the predictions, and the lower and upper
  limits of the cofidence intervals when `type = "subject_specific"`.

- nsim:

  numeric scalar giving the number of times to simulate the response
  variable.

- new_RE:

  logical; if `TRUE`, new random effects will be simulated, and new
  outcome data will be simulated by `simulate()` using these new random
  effect. Otherwise, the empirical Bayes estimates of the random effects
  from the fitted model will be used.

- acount_MLEs_var:

  logical; if `TRUE` it accounts for the variability of the maximum
  likelihood estimates (MLEs) by simulating a new value for the
  parameters from a multivariate normal distribution with mean the MLEs
  and covariance matrix the covariance matrix of the MLEs.

- sim_fun:

  a function based on which to simulate the response variable. This is
  relevant for non-standard models. The `simulate()` function also tries
  to extract this function from the `family` component of `object`. The
  function should have the following four arguments: `n` a numeric
  scalar denoting the number of observations to simulate, `mu` a numeric
  vector of means, `phis` a numeric vector of extra dispersion/scale
  parameters, and `eta_zi` a numeric vector for the zero-part of the
  model, if this is relevant.

- cores:

  the number of cores to use in the computation.

- ...:

  further arguments; currently none is used.

## Details

In generic terms, we assume that the mean of the outcome \\y_i\\ (\\i =
1, ..., n\\ denotes the subjects) conditional on the random effects is
given by the equation: \$\$g{E(y_i \| b_i)} = \eta_i = X_i \beta + Z_i
b_i,\$\$ where \\g(.)\\ denotes the link function, \\b_i\\ the vector of
random effects, \\\beta\\ the vector of fixed effects, and \\X_i\\ and
\\Z_i\\ the design matrices for the fixed and random effects,
respectively.

Argument `type_pred` of `predict()` specifies whether predictions will
be calculated in the link / linear predictor scale, i.e., \\\eta_i\\ or
in the response scale, i.e., \\g{E(y_i \| b_i)}\\.

When `type = "mean_subject"`, predictions are calculated using only the
fixed effects, i.e., the \\X_i \beta\\ part, where \\X_i\\ is evaluated
in `newdata`. This corresponds to predictions for the 'mean' subjects,
i.e., subjects who have random effects value equal to 0. Note, that in
the case of nonlinear link functions this does not correspond to the
averaged over the population predictions (i.e., marginal predictions).

When `type = "marginal"`, predictions are calculated using only the
fixed effects, i.e., the \\X_i \beta\\ part, where \\X_i\\ is evaluated
in `newdata`, but with \\\beta\\ coefficients the marginalized
coefficients obtain from
[`marginal_coefs`](https://drizopoulos.github.io/GLMMadaptive/reference/marginal_coefs.md).
These predictions will be marginal / population averaged predictions.

When `type = "zero_part"`, predictions are calculated for the logistic
regression of the extra zero-part of the model (i.e., applicable for
zero-inflated and hurdle models).

When `type = "subject_specific"`, predictions are calculated using both
the fixed- and random-effects parts, i.e., \\X_i \beta + Z_i b_i\\,
where \\X_i\\ and \\Z_i\\ are evaluated in `newdata`. Estimates for the
random effects of each subject are obtained as modes from the posterior
distribution \\\[b_i \| y_i; \theta\]\\ evaluated in `newdata` and with
\\theta\\ (denoting the parameters of the model, fixed effects and
variance components) replaced by their maximum likelihood estimates.

**Notes:** (i) When `se.fit = TRUE` and `type_pred = "response"`, the
standard errors returned are on the linear predictor scale, not the
response scale. (ii) When `se.fit = TRUE` and the model contains an
extra zero-part, no standard errors are computed when
`type = "mean_subject"`. (iii) When the model contains an extra
zero-part, `type = "marginal"` predictions are not yet implemented.

When `se.fit = TRUE` and `type = "subject_specific"`, standard errors
and confidence intervals for the subject-specific predictions are
obtained by a Monte Carlo scheme entailing three steps repeated `M`
times, namely

- Step I:

  Account for the variability of maximum likelihood estimates (MLES) by
  simulating a new value \\\theta^\*\\ for the parameters \\\theta\\
  from a multivariate normal distribution with mean the MLEs and
  covariance matrix the covariance matrix of the MLEs.

- Step II:

  Account for the variability in the random effects estimates by
  simulating a new value \\b_i^\*\\ for the random effects \\b_i\\ from
  the posterior distribution \\\[b_i \| y_i; \theta^\*\]\\. Because the
  posterior distribution does not have a closed-form, a
  Metropolis-Hastings algorithm is used to sample the new value
  \\b_i^\*\\ using as proposal distribution a multivariate Student's-t
  distribution with degrees of freedom `df`, centered at the mode of the
  posterior distribution \\\[b_i \| y_i; \theta\]\\ with \\\theta\\ the
  MLEs, and scale matrix the inverse Hessian matrix of the log density
  of \\\[b_i \| y_i; \theta\]\\ evaluated at the modes, but multiplied
  by `scale`. The `scale` and `df` parameters can be used to adjust the
  acceptance rate.

- Step III:

  The predictions are calculated using \\X_i \beta^\* + Z_i b_i^\*\\.

Argument `newdata2` can be used to calculate dynamic subject-specific
predictions. I.e., using the observed responses \\y_i\\ in `newdata`,
estimates of the random effects of each subject are obtained. For the
same subjects we want to obtain predictions in new covariates settings
for which no response data are yet available. For example, in a
longitudinal study, for a subject we have responses up to a follow-up
\\t\\ (`newdata`) and we want the prediction at \\t + \Delta t\\
(`newdata2`).

## Value

The estimated fixed and random effects, coefficients (this is similar as
in package **nlme**), confidence intervals fitted values (on the scale
on the response) and residuals.

## Author

Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>

## See also

[`mixed_model`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.md),
[`marginal_coefs`](https://drizopoulos.github.io/GLMMadaptive/reference/marginal_coefs.md)

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

fm1 <- mixed_model(fixed = y ~ year + group, random = ~ year | id, data = DF,
                   family = binomial())

head(coef(fm1))
#>   (Intercept)        year    group1
#> 1   -2.305583  0.54222964 -0.265744
#> 2   -1.905620 -0.12525835 -0.265744
#> 3   -2.395862 -0.14287663 -0.265744
#> 4   -2.379605  0.42518796 -0.265744
#> 5   -2.348230  0.01936688 -0.265744
#> 6   -2.056225  0.69488984 -0.265744
fixef(fm1)
#> (Intercept)        year      group1 
#>  -2.2391690   0.1968593  -0.2657440 
head(ranef(fm1))
#>   (Intercept)       year
#> 1 -0.06641381  0.3453703
#> 2  0.33354914 -0.3221176
#> 3 -0.15669315 -0.3397359
#> 4 -0.14043607  0.2283287
#> 5 -0.10906060 -0.1774924
#> 6  0.18294398  0.4980305


confint(fm1)
#>                  2.5 %   Estimate      97.5 %
#> (Intercept) -2.4728845 -2.2391690 -2.00545354
#> year         0.1529535  0.1968593  0.24076514
#> group1      -0.5673166 -0.2657440  0.03582867
confint(fm1, "var-cov")
#>                       2.5 %    Estimate     97.5 %
#> var.(Intercept)  0.22607344  0.50791169 1.14110835
#> cov.(Int)_year  -0.08018454 -0.03554827 0.07358216
#> var.year         0.15508238  0.16033127 0.20147617

head(fitted(fm1, "subject_specific"))
#>          1          2          3          4          5          6 
#> 0.09066165 0.15603499 0.83100175 0.96225413 0.97895242 0.97984947 
head(residuals(fm1, "marginal"))
#>          1          2          3          4          5          6 
#> -0.1805435 -0.1980669  0.6881459  0.6194025  0.5928527  0.5908654 

fm2 <- mixed_model(fixed = y ~ year * group, random = ~ year | id, data = DF,
                   family = binomial())

# likelihood ratio test between fm1 and fm2
anova(fm1, fm2)
#> 
#>         AIC     BIC  log.Lik  LRT df p.value
#> fm1 5203.04 5228.33 -2595.52                
#> fm2 5201.62 5231.12 -2593.81 3.43  1  0.0642
#> 

# the same but with a Wald test
anova(fm2, L = rbind(c(0, 0, 0, 1)))
#> 
#> Marginal Wald Tests Table
#> 
#> User-defined contrasts matrix:
#>  (Intr) year group1 yr:gr1
#>       0    0      0      1
#> 
#>    Chisq df Pr(>|Chi|)
#>   3.7927  1     0.0515
#> 
# }
```
