# GLMMadaptive Basics

## Generalized Linear Mixed Models - Theory

### Model Specification

Package **GLMMadaptive** provides a suit of functions for fitting and
post-processing mixed effects models for grouped/clustered outcomes
which have a distribution other than a normal distribution. In
particular, let y_i denote a vector of grouped/clustered outcome for the
i-th sample unit (i = 1, \ldots, n). The conditional distribution of y_i
given a vector of random effects b_i is assumed to be a member of the
extended exponential family, with linear predictor given by g\\E(y_i
\mid b_i)\\ = X_i \beta + Z_i b_i, where g(\cdot) denotes a monotonic
link function, X_i a design matrix for the fixed effects coefficients
\beta, and Z_i a design matrix for the random effects coefficients b_i.
Typically, matrix Z_i is assumed to be a subset of X_i. The random
effects are assumed to follow a normal distribution with mean 0 and
variance-covariance matrix D. In addition, the distribution \[y_i \mid
b_i\] may potentially have extra dispersion/shape parameters \phi.

### Estimation

The package focuses in settings in which the distribution \[y_i \mid
b_i\] is not normal and/or the link function g(\cdot) is not the
identity. In these settings, the estimation of the model is complicated
by the fact that the marginal log-likelihood function of the observed
y_i cannot be derived analytically. In particular, the log-likelihood
function has the form: \begin{array}{rcl} \ell(\theta) & = & \sum\_{i =
1}^n \log p(y_i; \theta)\\&&\\ & = & \sum\_{i = 1}^n \log \displaystyle
\int p(y_i \mid b_i; \theta) \\ p(b_i; \theta) \\ db_i, \end{array}
where \theta denotes the full parameter vector including the fixed
effects, the extra potential dispersion/shape parameters \phi and the
unique element of the covariance matrix D, and p(\cdot) denotes a
probability density or probability mass function. The integral in the
definition of \ell(\theta) does not have a closed-form solution, and
numerical approximations are required to obtain the maximum likelihood
estimates.

In the literature several approaches have been proposed to approximate
such integrals, and a nice overview is given in [Pinheiro and Chao
(2006)](https://doi.org/10.1198/106186006X96962). A typical approach to
approximate these integrals is the Laplace approximation. However, the
general consensus has been that in the standard but difficult cases of
binary/dichotomous data and count data with small counts and few
repeated measurements, the accuracy of this approximation is rather low.
Due to this fact, the general consensus is that the gold standard
numerical approximation method is the adaptive Gaussian quadrature rule
(*note: we focus here on maximum likelihood estimation; under the
Bayesian paradigm, approaches, such as, MCMC and Hamiltonian Monte Carlo
also provide accurate evaluation of the integrals*). This is more
computationally intensive but also more accurate. This package provides
an efficient implementation of the adaptive Gaussian quadrature rule,
allowing for multiple correlated random effects (e.g., random
intercepts, linear and quadratic random slopes) but currently a single
grouping factor (i.e., no nested or crossed random effects designs).

A hybrid optimization procedure is implemented starting with an EM
algorithm, treating the random effects as ‘missing data’, followed by a
direct optimization procedure with a quasi-Newton algorithm. Fine
control of this procedure is allowed with a series of control arguments.

## Generalized Linear Mixed Models - Practice

### Mixed Effects Logistic Regression

We illustrate the use of the package in the standard case of a mixed
effects logistic regression. That is, the distribution of \[y_i \mid
b_i\] is binomial, and the distribution of \[b_i\] multivariate normal.

We start by simulating some data for a binary longitudinal outcome:

``` r

set.seed(1234)
n <- 100 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 15 # maximum follow-up time

# we construct a data frame with the design: 
# everyone has a baseline measurement, and then measurements at random follow-up times
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ time, data = DF)

betas <- c(-2.13, -0.25, 0.24, -0.05) # fixed effects coefficients
D11 <- 0.48 # variance of random intercepts
D22 <- 0.1 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor
eta_y <- drop(X %*% betas + rowSums(Z * b[DF$id, ]))
# we simulate binary longitudinal data
DF$y <- rbinom(n * K, 1, plogis(eta_y))
```

We fit a mixed effects logistic regression for `y`, assuming random
intercepts for the random-effects part. The basic model-fitting function
in **GLMMadaptive** is called
[`mixed_model()`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.md),
and has four required arguments, namely `fixed` a formula for the fixed
effects, `random` a formula for the random effects, `family` a family
object specifying the type of response variable, and `data` a data frame
containing the variables in the previously mentioned formulas. Hence,
the call to fit the random intercepts logistic regression is:

``` r

fm1 <- mixed_model(fixed = y ~ sex * time, random = ~ 1 | id, data = DF,
                   family = binomial())
```

The summary method gives a detailed output of the model:

``` r

summary(fm1)
#> 
#> Call:
#> mixed_model(fixed = y ~ sex * time, random = ~1 | id, data = DF, 
#>     family = binomial())
#> 
#> Data Descriptives:
#> Number of Observations: 800
#> Number of Groups: 100 
#> 
#> Model:
#>  family: binomial
#>  link: logit 
#> 
#> Fit statistics:
#>    log.Lik      AIC      BIC
#>  -374.5197 759.0394 772.0653
#> 
#> Random effects covariance matrix:
#>               StdDev
#> (Intercept) 2.064622
#> 
#> Fixed effects:
#>                Estimate Std.Err z-value  p-value
#> (Intercept)     -2.6182  0.4448 -5.8862  < 1e-04
#> sexfemale       -1.3242  0.6436 -2.0575 0.039639
#> time             0.2851  0.0396  7.2089  < 1e-04
#> sexfemale:time   0.0402  0.0545  0.7372 0.461005
#> 
#> Integration:
#> method: adaptive Gauss-Hermite quadrature rule
#> quadrature points: 11
#> 
#> Optimization:
#> method: hybrid EM and quasi-Newton
#> converged: TRUE
```

We continue by checking the impact of the chosen number of quadrature
points to the parameters estimates and the log-likelihood value at
convergence. First, we refit the model with an increasing number of
quadrature points. The default when the number of random effects is
smaller or equal to two is 11 points. We fit then with 15, and 21
points:

``` r

fm1_q11 <- fm1
fm1_q15 <- update(fm1_q11, nAGQ = 15)
fm1_q21 <- update(fm1_q11, nAGQ = 21)

models <- list("nAGQ=11" = fm1_q11, "nAGQ=15" = fm1_q15, "nAGQ=21" = fm1_q21)
```

We now extract from the model the estimated parameter for the fixed
effects (using function
[`fixef()`](https://drizopoulos.github.io/GLMMadaptive/reference/methods.md)),
for the random effects, and the log-likelihood (using function
[`logLik()`](https://rdrr.io/r/stats/logLik.html)):

``` r

extract <- function (obj) {
    c(fixef(obj), "var_(Intercept)" = obj$D[1, 1], "logLik" = logLik(obj))
}

sapply(models, extract)
#>                       nAGQ=11       nAGQ=15      nAGQ=21
#> (Intercept)       -2.61816557   -2.61679077   -2.6168440
#> sexfemale         -1.32418678   -1.32349829   -1.3235049
#> time               0.28514773    0.28501148    0.2850151
#> sexfemale:time     0.04015372    0.04020003    0.0401995
#> var_(Intercept)    4.26266211    4.26351799    4.2636311
#> logLik          -374.51972108 -374.51983340 -374.5198046
```

We observe a rather stable model with virtually no differences between
the different choices of quadrature points.

We first compare the model with a simple logistic regression that does
not include any random effects:

``` r

km <- glm(y ~ sex * time, data = DF, family = binomial())

anova(fm1, km)
#> 
#>        AIC    BIC log.Lik    LRT df p.value
#> fm1 759.04 772.07 -374.52                  
#> km  919.11 937.85 -455.55 162.07  1 <0.0001
```

We obtain a highly significant p-value suggesting that there are
correlations in the data that cannot be ignored. *Note:* the
[`anova()`](https://drizopoulos.github.io/GLMMadaptive/reference/methods.md)
method that performs the likelihood ratio test calculates the p-value
using the standard \chi^2 distribution, here with one degree of freedom.
However, because the null hypothesis for testing variance parameters is
on the boundary of the corresponding parameter space, it would be more
appropriate to use a mixture of \chi^2 distributions.

We extend model `fm1` by also including a random slopes term; however,
we assume that the covariance between the random intercepts and random
slopes is zero. This is achieved by using the `||` symbol in the
specification of the `random` argument, i.e.,

``` r

fm2 <- mixed_model(fixed = y ~ sex * time, random = ~ time || id, data = DF,
                   family = binomial())
```

The likelihood ratio test between the two models is computed with
function
[`anova()`](https://drizopoulos.github.io/GLMMadaptive/reference/methods.md).
When two `"MixMod"` objects are provided, the function assumes that the
first object represents the model under the null hypothesis, and the
second object the model under the alternative, i.e.,

``` r

anova(fm1, fm2)
#> 
#>        AIC    BIC log.Lik  LRT df p.value
#> fm1 759.04 772.07 -374.52                
#> fm2 730.94 746.57 -359.47 30.1  1 <0.0001
```

The results suggest that we need the random slopes term. We continue by
testing whether the covariance between the random effects terms is zero.
The model under the alternative hypothesis is (results not shown):

``` r

fm3 <- mixed_model(fixed = y ~ sex * time, random = ~ time | id, data = DF,
                   family = binomial())
```

And again the likelihood ratio test is performed by a call to
[`anova()`](https://drizopoulos.github.io/GLMMadaptive/reference/methods.md)
(results not shown):

``` r

anova(fm2, fm3)
#> 
#>        AIC    BIC log.Lik  LRT df p.value
#> fm2 730.94 746.57 -359.47                
#> fm3 731.62 749.86 -358.81 1.32  1  0.2499
```

The results now suggest that indeed the covariance between the two
random effects terms is not statistically different from zero.

### Penalized Mixed Effects Poisson Regression

We continue our illustration with a Poisson mixed effects model. We
start again by simulating some data for a Poisson longitudinal outcome:

``` r

set.seed(1234)
n <- 100 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 15 # maximum follow-up time

# we construct a data frame with the design: 
# everyone has a baseline measurement, and then measurements at random follow-up times
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)

betas <- c(2.13, -0.25, 0.24, -0.05) # fixed effects coefficients
D11 <- 0.48 # variance of random intercepts

# we simulate random effects
b <- rnorm(n, sd = sqrt(D11))
# linear predictor
eta_y <- drop(X %*% betas + b[DF$id])
# we simulate Poisson longitudinal data
DF$y <- rpois(n * K, exp(eta_y))
```

We fit the mixed effects Poisson regression for `y` assuming random
intercepts for the random-effects part.

``` r

gm1 <- mixed_model(fixed = y ~ sex * time, random = ~ 1 | id, data = DF,
                   family = poisson())
```

As previously, the summary method gives a detailed output of the model:

``` r

summary(gm1)
#> 
#> Call:
#> mixed_model(fixed = y ~ sex * time, random = ~1 | id, data = DF, 
#>     family = poisson())
#> 
#> Data Descriptives:
#> Number of Observations: 800
#> Number of Groups: 100 
#> 
#> Model:
#>  family: poisson
#>  link: log 
#> 
#> Fit statistics:
#>    log.Lik      AIC      BIC
#>  -2782.306 5574.611 5587.637
#> 
#> Random effects covariance matrix:
#>                StdDev
#> (Intercept) 0.9781394
#> 
#> Fixed effects:
#>                Estimate Std.Err  z-value   p-value
#> (Intercept)      2.9447  0.3082   9.5545   < 1e-04
#> sexfemale       -0.8611  0.2699  -3.1906 0.0014199
#> time             0.2404  0.0016 151.3464   < 1e-04
#> sexfemale:time  -0.0511  0.0028 -18.3982   < 1e-04
#> 
#> Integration:
#> method: adaptive Gauss-Hermite quadrature rule
#> quadrature points: 11
#> 
#> Optimization:
#> method: hybrid EM and quasi-Newton
#> converged: TRUE
```

In settings with few subjects, repeated measurements of separation
effects, the package also allows to place a penalty in the fixed effects
regression coefficients \beta. The penalty/prior is in the form of a
Student’s t distribution with mean 0, scale parameter 1, and 3 degrees
of freedom, and it is placed in all \beta coefficients except from the
intercept. The penalized model can be fitted by setting argument
`penalized` to `TRUE`, i.e.,

``` r

gm2 <- mixed_model(fixed = y ~ sex * time, random = ~ 1 | id, data = DF,
                   family = poisson(), penalized = TRUE)
```

``` r

cbind('unpenalized' = fixef(gm1), 'penalized' = fixef(gm2))
#>                unpenalized   penalized
#> (Intercept)     2.94471505  2.94137865
#> sexfemale      -0.86111217 -0.85519416
#> time            0.24039444  0.24045750
#> sexfemale:time -0.05105926 -0.05114862
```

In this example we observe small differences between the penalized and
unpenalized models. The users have the option to alter the specification
of the Student’s t penalty by directly specifying the mean, scale and
degrees of freedom arguments of the distribution. For example, a ridge
penalty could be placed by setting the degrees of freedom to a high
value. The call in this case should be:

``` r

gm3 <- mixed_model(fixed = y ~ sex * time, random = ~ 1 | id, data = DF,
                   family = poisson(), 
                   penalized = list(pen_mu = 0, pen_sigma = 1, pen_df = 200))
```
