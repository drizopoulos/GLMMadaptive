# Family function for Negative Binomial Mixed Models

Specifies the information required to fit a Negative Binomial
generalized linear mixed model, using
[`mixed_model()`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.md).

## Usage

``` r
negative.binomial()
```

## Note

Currently only the log-link is implemented.

## Author

Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>

## Examples

``` r
# \donttest{
# simulate some data
set.seed(102)
dd <- expand.grid(f1 = factor(1:3), f2 = LETTERS[1:2], g = 1:30, rep = 1:15,
                  KEEP.OUT.ATTRS = FALSE)
mu <- 5 * (-4 + with(dd, as.integer(f1) + 4 * as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

gm1 <-  mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                    family = negative.binomial())

summary(gm1)
#> 
#> Call:
#> mixed_model(fixed = y ~ f1 * f2, random = ~1 | g, data = dd, 
#>     family = negative.binomial())
#> 
#> Data Descriptives:
#> Number of Observations: 2700
#> Number of Groups: 30 
#> 
#> Model:
#>  family: negative binomial
#>  link: log 
#> 
#> Fit statistics:
#>    log.Lik      AIC      BIC
#>  -10031.33 20078.66 20089.87
#> 
#> Random effects covariance matrix:
#>                 StdDev
#> (Intercept) 0.05344876
#> 
#> Fixed effects:
#>             Estimate Std.Err z-value   p-value
#> (Intercept)   1.6813  0.0719 23.3977   < 1e-04
#> f12           0.6208  0.0998  6.2215   < 1e-04
#> f13           0.9999  0.0994 10.0578   < 1e-04
#> f2B           1.6144  0.0990 16.3021   < 1e-04
#> f12:f2B      -0.4724  0.1394 -3.3900 0.0006989
#> f13:f2B      -0.7327  0.1390 -5.2696   < 1e-04
#> 
#> log(dispersion) parameter:
#>   Estimate Std.Err
#>    -0.7374  0.0277
#> 
#> Integration:
#> method: adaptive Gauss-Hermite quadrature rule
#> quadrature points: 11
#> 
#> Optimization:
#> method: hybrid EM and quasi-Newton
#> converged: TRUE 

# We do a likelihood ratio test with the Poisson mixed model
gm0 <- mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                    family = poisson())
#> Warning: Hessian matrix at convergence is not positive definite; unstable solution.
                    
anova(gm0, gm1)
#> 
#>          AIC      BIC   log.Lik      LRT df p.value
#> gm0 93397.20 93407.01 -46691.60                    
#> gm1 20078.66 20089.87 -10031.33 73320.54  1 <0.0001
#> 

# Define a cutom-made family function to be used with mixed_model()
# the required components are 'family', 'link', 'linkfun', 'linkinv' and 'log_dens';
# the extra functions 'score_eta_fun' and 'score_phis_fun' can be skipped and will 
# internally approximated using numeric derivatives (though it is better that you provide
# them).
my_negBinom <- function (link = "log") {
    stats <- make.link(link)
    log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
        # the log density function
        phis <- exp(phis)
        mu <- mu_fun(eta)
        log_mu_phis <- log(mu + phis)
        comp1 <- lgamma(y + phis) - lgamma(phis) - lgamma(y + 1)
        comp2 <- phis * log(phis) - phis * log_mu_phis
        comp3 <- y * log(mu) - y * log_mu_phis
        out <- comp1 + comp2 + comp3
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis, eta_zi) {
        # the derivative of the log density w.r.t. mu
        phis <- exp(phis)
        mu_phis <- mu + phis
        comp2 <- - phis / mu_phis
        comp3 <- y / mu - y / mu_phis
        # the derivative of mu w.r.t. eta (this depends on the chosen link function)
        mu.eta <- mu
        (comp2 + comp3) * mu.eta
    }
    score_phis_fun <- function (y, mu, phis, eta_zi) {
        # the derivative of the log density w.r.t. phis
        phis <- exp(phis)
        mu_phis <- mu + phis
        comp1 <- digamma(y + phis) - digamma(phis)
        comp2 <- log(phis) + 1 - log(mu_phis) - phis / mu_phis
        comp3 <- - y / mu_phis
        comp1 + comp2 + comp3
    }
    structure(list(family = "user Neg Binom", link = stats$name, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun,
                   score_phis_fun = score_phis_fun),
              class = "family")
}

fm <-  mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                   family = my_negBinom(), n_phis = 1, 
                   initial_values = list("betas" = poisson()))

summary(fm)
#> 
#> Call:
#> mixed_model(fixed = y ~ f1 * f2, random = ~1 | g, data = dd, 
#>     family = my_negBinom(), n_phis = 1, initial_values = list(betas = poisson()))
#> 
#> Data Descriptives:
#> Number of Observations: 2700
#> Number of Groups: 30 
#> 
#> Model:
#>  family: user Neg Binom
#>  link: log 
#> 
#> Fit statistics:
#>    log.Lik      AIC      BIC
#>  -10031.31 20078.61 20089.82
#> 
#> Random effects covariance matrix:
#>                 StdDev
#> (Intercept) 0.05215413
#> 
#> Fixed effects:
#>             Estimate Std.Err z-value   p-value
#> (Intercept)   1.6812  0.0718 23.4084   < 1e-04
#> f12           0.6208  0.0998  6.2225   < 1e-04
#> f13           0.9999  0.0994 10.0585   < 1e-04
#> f2B           1.6144  0.0990 16.3025   < 1e-04
#> f12:f2B      -0.4725  0.1394 -3.3905 0.0006977
#> f13:f2B      -0.7327  0.1390 -5.2696   < 1e-04
#> 
#> phi parameters:
#>       Estimate Std.Err
#> phi_1  -0.7374  0.0192
#> 
#> Integration:
#> method: adaptive Gauss-Hermite quadrature rule
#> quadrature points: 11
#> 
#> Optimization:
#> method: hybrid EM and quasi-Newton
#> converged: TRUE 
# }
```
