GLMMadaptive: Generalized Linear Mixed Models using Adaptive Gaussian Quadrature
================

<img src="man/figures/logo.png" height="205" align="right"/>

[![Travis-CI Build Status](https://travis-ci.org/drizopoulos/GLMMadaptive.svg?branch=master)](https://travis-ci.org/drizopoulos/GLMMadaptive) [![CRAN status](http://www.r-pkg.org/badges/version/GLMMadaptive)](https://cran.r-project.org/package=GLMMadaptive) [![](https://cranlogs.r-pkg.org/badges/grand-total/GLMMadaptive)](https://CRAN.R-project.org/package=GLMMadaptive) [![Rdoc](http://www.rdocumentation.org/badges/version/GLMMadaptive)](http://www.rdocumentation.org/packages/GLMMadaptive) [![Download counter](http://cranlogs.r-pkg.org/badges/GLMMadaptive)](https://cran.r-project.org/package=GLMMadaptive) 

Description
------------

This repository contains the source files for the R package <strong>GLMMadaptive</strong>. 
This package fits mixed effects models for grouped / repeated measurements data for which
the integral over the random effects in the definition of the marginal likelihood cannot
be solved analytically. The package approximates these integrals using the adaptive 
Gauss-Hermite quadrature rule.

Multiple random effects terms can be included for the grouping factor (e.g., random 
intercepts, random linear slopes, random quadratic slopes), but currently only a single
grouping factor is allowed.

Basic Features
------------

- The package contains a single model-fitting function named `mixed_model()` with four 
required arguments, `fixed` a formula for the fixed effects, `random` a formula for the
random effects, `family` a family object specifying the type of response variable, and 
`data` a data frame containing the variables in the previously mentioned formulas.

- Methods for standard generics are provided, i.e., `coef()`, `fixef()`, `ranef()`, 
`vcov()`, `logLik()`, `summary()`, `anova()`, `confint()`, `fitted()`, `residuals()`, 
`predict()`, and `simulate()`.

- Negative binomial mixed models can be fitted using the `negative.binomial()` family 
object.

- Zero-inflated Poisson and negative binomial models using the `zi.poisson()` and 
`zi.negative.binomial()` family objects.

- Users may also specify their own log-density function for the repeated measurements 
response variable, and the internal algorithms will take care of the optimization.

- Calculates the marginalized coefficients using the idea of Hedeker et al. (2017) using 
function `marginal_coefs()`.

- Predictions with confidence interval for constructing effects plots are provided by 
function `effectPlotData()`.

Basic Use
------------

Let `y` denote a grouped/clustered outcome, `g` denote the grouping factor, and `x1` and
`x2` covariates. A mixed effects model with `y` as outcome, `x1` and `x2` as fixed effect,
and random intercepts is fitted with the code:
```r
fm <- mixed_model(y ~ x1 + x2, random = ~ 1 | g, data = DF)

summary(fm)
```

In the `data` argument we provide the data frame `DF`, which contains the aforementioned 
variables. To include in the random-effects part intercepts and `x1`, we update the call
to `mixed_model()` as
```r
gm <- mixed_model(y ~ x1 + x2, random = ~ x1 | g, data = DF)

summary(gm)
```


 
Installation
------------

The development version of the package can be installed from GitHub using the **devtools**
package:
```r
devtools::install_github("drizopoulos/GLMMadaptive")
```

and with vignettes
```r
devtools::install_github("drizopoulos/GLMMadaptive", build_vignettes = TRUE)
```

Hex-sticker courtesy of Greg Papageorgiou [@gr_papageorgiou](https://twitter.com/gr_papageorgiou).