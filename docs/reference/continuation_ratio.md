# Functions to Set-Up Data for a Continuation Ratio Mixed Model

Data set-up and calculation of marginal probabilities from a
continuation ratio model

## Usage

``` r
cr_setup(y, direction = c("forward", "backward"))

cr_marg_probs(eta, direction = c("forward", "backward"))
```

## Arguments

- y:

  a numeric vector denoting the ordinal response variable.

- direction:

  character string specifying the direction of the continuation ratio
  model; `"forward"` corresponds to a discrete hazard function.

- eta:

  a numeric matrix of the linear predictor, with columns corresponding
  to the different levels of the ordinal response.

## Author

Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>

Frank Harrell

## Note

Function `cr_setup()` is based on the `cr.setup()` function from package
**rms**.

## Examples

``` r
n <- 300 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 15 # maximum follow-up time

# we constuct a data frame with the design: 
# everyone has a baseline measurment, and then measurements at random follow-up times
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)[, -1]
Z <- model.matrix(~ 1, data = DF)

thrs <- c(-1.5, 0, 0.9) # thresholds for the different ordinal categories
betas <- c(-0.25, 0.24, -0.05) # fixed effects coefficients
D11 <- 0.48 # variance of random intercepts
D22 <- 0.1 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))[, 1, drop = FALSE]
# linear predictor
eta_y <- drop(X %*% betas + rowSums(Z * b[DF$id, , drop = FALSE]))
# linear predictor for each category
eta_y <- outer(eta_y, thrs, "+")
# marginal probabilities per category
mprobs <- cr_marg_probs(eta_y)
# we simulate ordinal longitudinal data
DF$y <- unname(apply(mprobs, 1, sample, x = ncol(mprobs), size = 1, replace = TRUE))

# If you want to simulate from the backward formulation of the CR model, you need to
# change `eta_y <- outer(eta_y, thrs, "+")` to `eta_y <- outer(eta_y, rev(thrs), "+")`,
# and `mprobs <- cr_marg_probs(eta_y)` to `mprobs <- cr_marg_probs(eta_y, "backward")`

#################################################

# prepare the data
# If you want to fit the CR model under the backward formulation, you need to change
# `cr_vals <- cr_setup(DF$y)` to `cr_vals <- cr_setup(DF$y, "backward")`
cr_vals <- cr_setup(DF$y)
cr_data <- DF[cr_vals$subs, ]
cr_data$y_new <- cr_vals$y
cr_data$cohort <- cr_vals$cohort

# fit the model
fm <- mixed_model(y_new ~ cohort + sex * time, random = ~ 1 | id, 
                  data = cr_data, family = binomial())

summary(fm)
#> 
#> Call:
#> mixed_model(fixed = y_new ~ cohort + sex * time, random = ~1 | 
#>     id, data = cr_data, family = binomial())
#> 
#> Data Descriptives:
#> Number of Observations: 4279
#> Number of Groups: 300 
#> 
#> Model:
#>  family: binomial
#>  link: logit 
#> 
#> Fit statistics:
#>    log.Lik      AIC      BIC
#>  -2490.064 4994.128 5020.055
#> 
#> Random effects covariance matrix:
#>                StdDev
#> (Intercept) 0.6983969
#> 
#> Fixed effects:
#>                Estimate Std.Err  z-value p-value
#> (Intercept)     -1.6528  0.1136 -14.5440 < 1e-04
#> cohorty>=2       1.4226  0.0878  16.1958 < 1e-04
#> cohorty>=3       2.4467  0.1320  18.5331 < 1e-04
#> sexfemale       -0.1741  0.1395  -1.2476 0.21218
#> time             0.2503  0.0132  19.0206 < 1e-04
#> sexfemale:time  -0.0661  0.0164  -4.0306 < 1e-04
#> 
#> Integration:
#> method: adaptive Gauss-Hermite quadrature rule
#> quadrature points: 11
#> 
#> Optimization:
#> method: EM
#> converged: TRUE 
```
