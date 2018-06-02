GLMMadaptive: Generalized Linear Mixed Models using Adaptive Gaussian Quadrature
================
[![Travis-CI Build Status](https://travis-ci.org/drizopoulos/GLMMadaptive.svg?branch=master)](https://travis-ci.org/drizopoulos/GLMMadaptive)

Description
------------

This repository contains the source files for the R package <strong>GLMMadaptive</strong>. 
This package fits mixed effects models for grouped / repeated measurements data with 
nonlinear link functions under the maximum likelihood approach. The integrals over the 
random effects are approximated with an adaptive Gaussian quadrature rule.

Multiple random effects terms can be included for the grouping factor (e.g., random 
intercepts, random linear slopes, random quadratice slopes), but currently only a single
grouping factor is allowed.

Basic Features
------------

- The package contains a sinle model-fitting function named `mixed_model()`.

- Methods for standard generics are provided, i.e., `coef()`, `fixef()`, `ranef()`, 
`vcov()`, `logLik()`, `summary()`, `anova()`, `confint()`, `fitted()`, `residuals()`, 
and `predict()`.

- Negative binomial mixed models can be fitted using the `negative.binomial()` family object.

- Users may also specify their own log-density function for the repeated measurements 
response variable, and the internal algorithms will take care of the optimization.

- Calculates the marginalized coefficients using the idea of Hedeker et al. (2017) using 
function `marginal_coefs()`.

- Predictions with confidence interval for constructing effects plots are provided by 
function `effectPlotData()`.
