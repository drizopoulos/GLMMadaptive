# Family functions for Student's-t, Beta, Zero-Inflated and Hurdle Poisson and Negative Binomial, Hurdle Log-Normal, Hurdle Beta, Gamma, and Censored Normal Mixed Models

Specifies the information required to fit a Beta, zero-inflated and
hurdle Poisson, zero-inflated and hurdle Negative Binomial, a hurdle
normal and a hurdle Beta mixed-effects model, using
[`mixed_model()`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.md).

## Usage

``` r
students.t(df = stop("'df' must be specified"), link = "identity")
beta.fam()
zi.poisson()
zi.binomial()
zi.negative.binomial()
hurdle.poisson()
hurdle.negative.binomial()
hurdle.lognormal()
hurdle.beta.fam()
unit.lindley()
beta.binomial(link = "logit")
Gamma.fam()
censored.normal()
```

## Arguments

- link:

  name of the link function.

- df:

  the degrees of freedom of the Student's t distribution.

## Note

Currently only the log-link is implemented for the Poisson, negative
binomial and Gamma models, the logit link for the beta and hurdle beta
models and the identity link for the log-normal model.

## Examples

``` r
# simulate some data from a negative binomial model
set.seed(102)
dd <- expand.grid(f1 = factor(1:3), f2 = LETTERS[1:2], g = 1:30, rep = 1:15,
                  KEEP.OUT.ATTRS = FALSE)
mu <- 5*(-4 + with(dd, as.integer(f1) + 4 * as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

# Fit a zero-inflated Poisson model, with only fixed effects in the 
# zero-inflated part
fm1 <- mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                  family = zi.poisson(), zi_fixed = ~ 1)

summary(fm1)
#> 
#> Call:
#> mixed_model(fixed = y ~ f1 * f2, random = ~1 | g, data = dd, 
#>     family = zi.poisson(), zi_fixed = ~1)
#> 
#> Data Descriptives:
#> Number of Observations: 2700
#> Number of Groups: 30 
#> 
#> Model:
#>  family: zero-inflated poisson
#>  link: log 
#> 
#> Fit statistics:
#>    log.Lik      AIC      BIC
#>  -38656.66 77329.32 77340.53
#> 
#> Random effects covariance matrix:
#>                StdDev
#> (Intercept) 0.1206881
#> 
#> Fixed effects:
#>             Estimate Std.Err  z-value p-value
#> (Intercept)   2.0369  0.0301  67.7765 < 1e-04
#> f12           0.5228  0.0253  20.6624 < 1e-04
#> f13           0.9018  0.0239  37.7840 < 1e-04
#> f2B           1.3992  0.0224  62.5730 < 1e-04
#> f12:f2B      -0.3802  0.0282 -13.4937 < 1e-04
#> f13:f2B      -0.6665  0.0268 -24.9152 < 1e-04
#> 
#> Zero-part coefficients:
#>             Estimate Std.Err  z-value p-value
#> (Intercept)  -1.4413  0.0489 -29.4463 < 1e-04
#> 
#> Integration:
#> method: adaptive Gauss-Hermite quadrature rule
#> quadrature points: 11
#> 
#> Optimization:
#> method: hybrid EM and quasi-Newton
#> converged: TRUE 

# \donttest{
# We extend the previous model allowing also for a random intercept in the
# zero-inflated part
fm2 <- mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                  family = zi.poisson(), zi_fixed = ~ 1, zi_random = ~ 1 | g)

# We do a likelihood ratio test between the two models
anova(fm1, fm2)
#> 
#>          AIC      BIC   log.Lik  LRT df p.value
#> fm1 77329.32 77340.53 -38656.66                
#> fm2 77331.95 77345.96 -38655.97 1.37  2  0.5047
#> 

#############################################################################
#############################################################################

# The same as above but with a negative binomial model
gm1 <- mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                  family = zi.negative.binomial(), zi_fixed = ~ 1)

summary(gm1)
#> 
#> Call:
#> mixed_model(fixed = y ~ f1 * f2, random = ~1 | g, data = dd, 
#>     family = zi.negative.binomial(), zi_fixed = ~1)
#> 
#> Data Descriptives:
#> Number of Observations: 2700
#> Number of Groups: 30 
#> 
#> Model:
#>  family: zero-inflated negative binomial
#>  link: log 
#> 
#> Fit statistics:
#>    log.Lik     AIC      BIC
#>  -10031.15 20080.3 20092.91
#> 
#> Random effects covariance matrix:
#>                 StdDev
#> (Intercept) 0.04718218
#> 
#> Fixed effects:
#>             Estimate Std.Err z-value    p-value
#> (Intercept)   1.6877  0.0738 22.8762    < 1e-04
#> f12           0.6210  0.0994  6.2470    < 1e-04
#> f13           1.0012  0.0992 10.0975    < 1e-04
#> f2B           1.6138  0.0987 16.3563    < 1e-04
#> f12:f2B      -0.4724  0.1388 -3.4028 0.00066693
#> f13:f2B      -0.7344  0.1386 -5.2979    < 1e-04
#> 
#> Zero-part coefficients:
#>             Estimate Std.Err z-value  p-value
#> (Intercept)  -5.0279   2.785 -1.8053 0.071024
#> 
#> log(dispersion) parameter:
#>   Estimate Std.Err
#>    -0.7205  0.0545
#> 
#> Integration:
#> method: adaptive Gauss-Hermite quadrature rule
#> quadrature points: 11
#> 
#> Optimization:
#> method: hybrid EM and quasi-Newton
#> converged: TRUE 

# We do a likelihood ratio test between the Poisson and negative binomial models
anova(fm1, gm1)
#> 
#>          AIC      BIC   log.Lik      LRT df p.value
#> fm1 77329.32 77340.53 -38656.66                    
#> gm1 20080.30 20092.91 -10031.15 57251.01  1 <0.0001
#> 
# }
```
