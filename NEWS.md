# GLMMadaptive 0.8.5

## General
* Added the `zi.binomial()` family object for fitting zero-inflated Binomial mixed-effects models.

# GLMMadaptive 0.8.0

## General
* Added the `beta.binomial()` family object for fitting Beta-Binomial mixed-effects models.

* Added the `Gamma.fam()` family object for fitting Gamma mixed-effects models.

* Added the `censored.normal()` family object for fitting linear mixed-effects models with right and left censored outcomes.

## Minor
* Corrected typos in the vignettes and the help files.


# GLMMadaptive 0.6.5

## General
* Added a `weights` argument in `mixed_model()`. These are simple multipliers on the log-likelihood contributions of each group/cluster, i.e., we presume that there are multiple replicates of each group/cluster denoted by the weights.

* The internal implementation of the `negative.binomial()` family is now considerably faster.

* More stable and accurate calculations using the **matrixStats** package.


# GLMMadaptive 0.6.0

## General
* Mixed models for ordinal clustered outcomes using the continuation ratio model, with forward and backward formulation.

* The new function `VIF()` calculates variance inflation factors for mixed models fitted in the package.

* The CRAN version of the package has now only three basic vignettes to conform with the time restrictions of CRAN checks. The whole list of available vignettes can be found in the website of the package [https://drizopoulos.github.io/GLMMadaptive/](https://drizopoulos.github.io/GLMMadaptive/).

## Minor
* Corrected a typo in `anova.MixMod()` that was reporting the AIC as the BIC of the fitted model.

* The offset is now passed in the calculation of initial values.

# GLMMadaptive 0.5.0

## General
* Added support for the **DHARMa** package.

* The new vignette *Goodness of Fit for MixMod Objects* describes how to check the fit 
of mixed models fitted by **GLMMadaptive**.

* Added support for the **effects** package.

* There is a new section in the vignette *Methods for MixMod Objects* illustrating the use of the **effects** package.

* Function `marginal_coefs()` has a faster implementation. Compared to the previous implementation the results will be slightly different.

* The optimizer `nlminb()` can now also be invoked using the new control argument `optimizer`; default is `"optim"` corresponding to function `optim()`.

* The new vignette *Optimization and Numerical Integration in GLMMadaptive* describes how to control the optimization and numerical integration procedures in the package.

# GLMMadaptive 0.4.0

## General
* The `predict()` method now works for zero-inflated and hurdle models.

* Hurdle Beta mixed effects models are now available using the `hurdle.beta.fam` family
object.

* The new function `scoring_rules()` calculates proper scoring rule for subject-specific
predictions from mixed models for categorical data.

* Added support for the **emmeans** package.

# GLMMadaptive 0.3.0

## General
* Hurdle Poisson and negative binomial models are now implemented using the family objects `hurdle.poisson` and `hurdle.negative.binomial`, respectively.

* added S3 methods for the `terms()`, `model.frame()` and `model.matrix()` generics in order to
work with the **multcomp** package.

* A new vignette illustrating multiple comparisons with the **multcomp** package.

* Methods `vcov()`, `summary()`, `confint()`, `anova()`, `marginal_coefs()`, `effectPlotData()`, `predict()`, and `simulate()` have gained the logical argument `sandwich` to invoke the use of robust/sandwich standard errors in the calculations.

# GLMMadaptive 0.2.0

## General
* Zero-inflated Poisson and negative binomial models are now implemented using the family objects `zi.poisson()` and `zi.negative.binomial()`, respectively. In addition, taking into advantage of the fact that users can specify their own log density functions for the outcome, two-part / hurdle model can also be implemented. 

* A new vignette illustrates how the zero-inflated models can be fitted.

* The `predict()` method is now fully available. It calculates predictions, and standard errors for models returned by `mixed_model()` at three levels:
     + "mean subject": only the fixed effects part corresponding to predictions for the average subject (but not population averaged predictions in case of nonlinear link functions).
     
     + "marginal": predictions using the marginalized coefficients that correspond to population averaged predictions.
     
     + "subject specific": predictions at the subject level. These can be also calculated for subjects not originally in the dataset (i.e., estimates of the random effects are internally obtained).

* The `simulate()` method is available to simulate data from fitted mixed models. This can be used for instance to perform replication / posterior predictive checks.

# GLMMadaptive 0.1.6

## General
* Added vignettes.

# GLMMadaptive 0.1.3

## General
* First version of the package.

